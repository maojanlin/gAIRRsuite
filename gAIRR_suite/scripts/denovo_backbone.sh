raw_seq_path=$1
contig_path=$2
contig_check_path=$3
flanking_result_path=$4
allele_name=$5
cluster_num=$6
path_SPAdes=$7
thread=$8


# get the script directory
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


mkdir -p ${contig_path} 
echo "[AIRRCall] [FLANKING SEQUENCE] Denovo assemble the backbones..."
if [ -f "${contig_path}spades_log.log" ] ; then
    rm ${contig_path}spades_log.log
fi
for ((cluster_id=0; cluster_id <${cluster_num}; cluster_id++ ))
do
    echo "[SPAdes] group ${cluster_id} assembled..."
    spades.py -1 ${raw_seq_path}${allele_name}_${cluster_id}_read_R1.fasta \
                           -2 ${raw_seq_path}${allele_name}_${cluster_id}_read_R2.fasta \
                           --only-assembler -t ${thread} \
                           -o ${contig_path}${allele_name}_${cluster_id} \
                           >> ${contig_path}spades_log.log
    python3 ${script_dir}/copy_1st_contig.py -fc ${contig_path}${allele_name}_${cluster_id}/contigs.fasta \
                                             -id ${cluster_id} \
                                             -map ${raw_seq_path}/${allele_name}_group_allele_map.txt \
                                             -fo ${contig_path}${allele_name}_contig_${cluster_id}.fasta
done

mkdir -p ${contig_check_path}
echo "[AIRRCall] [FLANKING SEQUENCE] Checking contigs..."
if [ -f "${contig_check_path}bwa_log.log" ] ; then
    rm ${contig_check_path}bwa_log.log
fi
if [ -f "${flanking_result_path}/assembly_call.txt" ] ; then
    rm ${flanking_result_path}/assembly_call.txt
fi
if [ -f "${flanking_result_path}/flanking_contigs.fasta" ] ; then
    rm ${flanking_result_path}/flanking_contigs.fasta
fi
if [ -f "${flanking_result_path}/flanking_size200.fasta" ] ; then
    rm ${flanking_result_path}/flanking_size200.fasta
fi
if [ -f "${flanking_result_path}/flank_region.txt" ] ; then
    rm ${flanking_result_path}/flank_region.txt
fi
for ((cluster_id=0; cluster_id <${cluster_num}; cluster_id++ ))
do
    echo "checking contig ${cluster_id}"
    bwa index     ${contig_path}/${allele_name}_contig_${cluster_id}.fasta \
                  2>> ${contig_check_path}bwa_log.log
    if [ ${allele_name} == "TCRD_plusHep" ] || [ ${allele_name} == "BCRD_plusHep" ]; then
        echo "[AIRRCall] Adjust BWA parameters for shorter alleles..."
        bwa mem -t ${thread} -T 10 ${contig_path}/${allele_name}_contig_${cluster_id}.fasta  \
                             ${raw_seq_path}/${allele_name}_${cluster_id}_allele.fasta \
                             >   ${contig_check_path}/align_${allele_name}_${cluster_id}.sam \
                             2>> ${contig_check_path}/bwa_log.log
    else
        bwa mem -t ${thread} ${contig_path}/${allele_name}_contig_${cluster_id}.fasta  \
                             ${raw_seq_path}/${allele_name}_${cluster_id}_allele.fasta \
                             >   ${contig_check_path}/align_${allele_name}_${cluster_id}.sam \
                             2>> ${contig_check_path}/bwa_log.log
    fi
    python3 ${script_dir}/parse_bwa_sam.py -fs ${contig_check_path}/align_${allele_name}_${cluster_id}.sam \
                                           -fc ${contig_path}/${allele_name}_contig_${cluster_id}.fasta \
                                           -td 0 \
                                           -fo ${flanking_result_path}/assembly_call.txt \
                                           -fr ${flanking_result_path}/flanking_contigs.fasta \
                                           --cluster_id ${cluster_id} \
                                           -fsize 200 \
                                           -frs ${flanking_result_path}/flanking_size200.fasta \
                                           >>   ${flanking_result_path}/flank_region.txt
done
