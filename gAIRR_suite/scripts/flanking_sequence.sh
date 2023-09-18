# input path parameter
outer_dir=$1/$4_$2_flanking/
allele_name=$2
allele_path=$3
read_path_1=$5
read_path_2=$6
path_SPAdes=$7
thread=$8

# get the script directory
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


# working sub directories
raw_seq_path=${outer_dir}group_allele_reads/
contig_path=${outer_dir}asm_contigs/
contig_check_path=${outer_dir}asm_check/
haplotype_sam_path=${outer_dir}haplotype_sam/
flanking_result_path=${outer_dir}flanking_result/

# building the environment and grouping the reads
mkdir -p ${flanking_result_path}
mkdir -p ${raw_seq_path}
echo "[AIRRCall] [FLANKING SEQUENCE] Grouping alleles and reads..."
cluster_num=$( python3 ${script_dir}/group_allele_reads.py -fp  $1/$4_$2/allele_support_reads.pickle \
                                                           -fa  ${allele_path} \
                                                           -fr1 ${read_path_1} \
                                                           -fr2 ${read_path_2} \
                                                           -name ${allele_name} \
                                                           -fod ${raw_seq_path} )

echo "[AIRRCall] [FLANKING SEQUENCE] Get the backbone flanking sequences..."
${script_dir}/denovo_backbone.sh ${raw_seq_path} ${contig_path} ${contig_check_path} ${flanking_result_path} ${allele_name} ${cluster_num} ${path_SPAdes} ${thread}

mkdir -p $haplotype_sam_path
echo "[AIRRCall] [FLANKING SEQUENCE] Align short reads to the backbones..."
bwa index ${flanking_result_path}/flanking_contigs.fasta
bwa mem -t ${thread} ${flanking_result_path}/flanking_contigs.fasta ${read_path_1} ${read_path_2} > ${haplotype_sam_path}/bwa_reads_to_flanking.sam

len_extend=200
if [ -f "${flanking_result_path}/flanking_haplotypes.fasta" ] ; then
    rm ${flanking_result_path}/flanking_haplotypes.fasta
fi
python3 ${script_dir}/shrink_sam_to_range.py -fs ${haplotype_sam_path}/bwa_reads_to_flanking.sam \
                                             -fc ${flanking_result_path}/flanking_contigs.fasta \
                                             -fr ${flanking_result_path}/flank_region.txt \
                                             -ext ${len_extend} \
                                             -foh ${flanking_result_path}/flanking_haplotypes.fasta \
                                             -foc ${haplotype_sam_path}/flanking_contigs_extend_${len_extend}.fasta \
                                             >    ${haplotype_sam_path}bwa_flanking_cover_link.rpt \
                                             2>   ${haplotype_sam_path}bwa_flanking_all.txt

echo "[AIRRCall] [FLANKING SEQUENCE] Finish flanking sequence calling."

