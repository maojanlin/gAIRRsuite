# input path parameter
outer_dir=$1/$4_$2_flanking/
allele_name=$2
allele_path=$3
read_path_1=$5
read_path_2=$6

# working sub directories
raw_seq_path=${outer_dir}group_allele_reads/
contig_path=${outer_dir}asm_contigs/
contig_check_path=${outer_dir}asm_check/
haplotype_sam_path=${outer_dir}haplotype_sam/
flanking_result_path=${outer_dir}flanking_result/

# building the environment and grouping the reads
mkdir -p ${raw_seq_path}
echo "[AIRRCall] [FLANKING SEQUENCE] Grouping alleles and reads..."
cluster_num=$( python3 group_allele_reads.py -fp  $1/$4_$2/allele_support_reads.pickle \
                                           -fa  ${allele_path} \
                                           -fr1 ${read_path_1} \
                                           -fr2 ${read_path_2} \
                                           -name ${allele_name} \
                                           -fod ${raw_seq_path} )

mkdir -p ${contig_path} 
rm ${contig_path}spades_log.log
for ((cluster_id=0; cluster_id <${cluster_num}; cluster_id++ ))
do
    echo "[SPAdes] group ${cluster_id} assembled..."
    python3 ../SPAdes-3.11.1-Darwin/bin/spades.py -1 ${raw_seq_path}${allele_name}_${cluster_id}_read_R1.fasta \
                                                  -2 ${raw_seq_path}${allele_name}_${cluster_id}_read_R2.fasta \
                                                  --only-assembler -t 8 \
                                                  -o ${contig_path}${allele_name}_${cluster_id} \
                                                  >> ${contig_path}spades_log.log
    python3 copy_1st_contig.py -fc ${contig_path}${allele_name}_${cluster_id}/contigs.fasta \
                               -fo ${contig_path}${allele_name}_contig_${cluster_id}.fasta
done

mkdir -p ${contig_check_path}
echo "[AIRRCall] [FLANKING SEQUENCE] Checking contigs..."
rm ${contig_check_path}bwa_log.log
for ((cluster_id=0; cluster_id <${cluster_num}; cluster_id++ ))
do
    bwa index     ${contig_path}/${allele_name}_contig_${cluster_id}.fasta \
                  2>> ${contig_check_path}bwa_log.log
    bwa mem -t 16 ${contig_path}/${allele_name}_contig_${cluster_id}.fasta  \
                  ${raw_seq_path}/${allele_name}_${cluster_id}_allele.fasta \
                  >   ${contig_check_path}/align_${allele_name}_${cluster_id}.sam \
                  2>> ${contig_check_path}/bwa_log.log
done

mkdir -p $haplotype_sam_path
mkdir -p $flanking_result_path
echo "[AIRRCall] [FLANKING SEQUENCE] Cropping and Haplotyping the assembly..."
rm ${haplotype_sam_path}/bwa_read_to_contig_log.log
rm ${flanking_result_path}/assembly_call.txt
rm ${flanking_result_path}/flanking_contigs.fasta
rm ${flanking_result_path}/flanking_size200.fasta
rm ${flanking_result_path}/corrected_contigs.fasta
for ((cluster_id=0; cluster_id <${cluster_num}; cluster_id++ ))
do
    echo "align reads to contig_${cluster_id}."
    flank_region=$( python3 parse_bwa_sam.py -fs ${contig_check_path}/align_${allele_name}_${cluster_id}.sam \
                                             -fc ${contig_path}/${allele_name}_contig_${cluster_id}.fasta \
                                             -td 0 \
                                             -fo ${flanking_result_path}/assembly_call.txt \
                                             -fr ${flanking_result_path}/flanking_contigs.fasta \
                                             --cluster_id ${cluster_id}-1 \
                                             -fsize 200 \
                                             -frs ${flanking_result_path}/flanking_size200.fasta)
    # if parse_bwa_sam.py report region, then do the 2nd round assembly
    if [ ! -z ${flank_region} ]  # if $flank_region exist
    then
        # align reads back to the SPAdes_assembly contigs as "realign.sam"
        bwa mem -t 16 ${contig_path}/${allele_name}_contig_${cluster_id}.fasta \
                      ${raw_seq_path}/${allele_name}_${cluster_id}_read_R1.fasta \
                      ${raw_seq_path}/${allele_name}_${cluster_id}_read_R2.fasta \
                      > ${haplotype_sam_path}/${allele_name}_realign_${cluster_id}.sam \
                      2>> ${haplotype_sam_path}/bwa_read_to_contig_log.log
        # parse the "realign.sam", and pop the perfect match to the 1st_round_contigs to generate the alternative remain_P1.fasta and remain_P2.fasta
        python3 parse_contig_realign.py -fs ${haplotype_sam_path}/${allele_name}_realign_${cluster_id}.sam \
                                        -fo ${haplotype_sam_path}/${allele_name}_remain_${cluster_id}.fasta \
                                        -it_region ${flank_region} \
                                        --contig_file ${contig_path}/${allele_name}_contig_${cluster_id}.fasta \
                                        --allele_file ${raw_seq_path}/${allele_name}_${cluster_id}_allele.fasta \
                                        --corrected_contig_output_file ${flanking_result_path}/corrected_contigs.fasta \
                                        > ${haplotype_sam_path}/${allele_name}_realign_${cluster_id}.rpt
    fi
done
echo "[AIRRCall] [FLANKING SEQUENCE] Finish flanking sequence calling."

