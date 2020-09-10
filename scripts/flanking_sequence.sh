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
mkdir -p $flanking_result_path
mkdir -p ${raw_seq_path}
echo "[AIRRCall] [FLANKING SEQUENCE] Grouping alleles and reads..."
cluster_num=$( python3 group_allele_reads.py -fp  $1/$4_$2/allele_support_reads.pickle \
                                           -fa  ${allele_path} \
                                           -fr1 ${read_path_1} \
                                           -fr2 ${read_path_2} \
                                           -name ${allele_name} \
                                           -fod ${raw_seq_path} )

echo "[AIRRCall] [FLANKING SEQUENCE] Get the backbone flanking sequences..."
./denovo_backbone.sh ${raw_seq_path} ${contig_path} ${contig_check_path} ${flanking_result_path} ${allele_name} ${cluster_num}

mkdir -p $haplotype_sam_path
echo "[AIRRCall] [FLANKING SEQUENCE] Align short reads to the backbones..."
bwa index ${flanking_result_path}/flanking_contigs.fasta
bwa mem -t 16 ${flanking_result_path}/flanking_contigs.fasta ${read_path_1} ${read_path_2} > ${haplotype_sam_path}/bwa_reads_to_flanking.sam

len_extend=200
rm ${flanking_result_path}/flanking_haplotypes.fasta
python3 shrink_sam_to_range.py -fs ${haplotype_sam_path}/bwa_reads_to_flanking.sam \
                               -fc ${flanking_result_path}/flanking_contigs.fasta \
                               -fr ${flanking_result_path}/flank_region.txt \
                               -ext ${len_extend} \
                               -foh ${flanking_result_path}/flanking_haplotypes.fasta \
                               -foc ${haplotype_sam_path}/flanking_contigs_extend_${len_extend}.fasta \
                               >    ${haplotype_sam_path}bwa_flanking_cover_link.rpt \
                               2>   ${haplotype_sam_path}bwa_flanking_all.txt

#mkdir -p $haplotype_sam_path
#echo "[AIRRCall] [FLANKING SEQUENCE] Cropping and Haplotyping the assembly..."
#rm ${haplotype_sam_path}/bwa_read_to_contig_log.log
#rm ${flanking_result_path}/corrected_contigs.fasta
#for ((cluster_id=0; cluster_id <${cluster_num}; cluster_id++ ))
#do
#    # if parse_bwa_sam.py report region, then do the 2nd round assembly
#    if [ ! -z ${flank_region} ]  # if $flank_region exist
#    then
#        echo "align reads to contig_${cluster_id}."
#        # align reads back to the SPAdes_assembly contigs as "realign.sam"
#        bwa mem -t 16 ${contig_path}/${allele_name}_contig_${cluster_id}.fasta \
#                      ${raw_seq_path}/${allele_name}_${cluster_id}_read_R1.fasta \
#                      ${raw_seq_path}/${allele_name}_${cluster_id}_read_R2.fasta \
#                      > ${haplotype_sam_path}/${allele_name}_realign_${cluster_id}.sam \
#                      2>> ${haplotype_sam_path}/bwa_read_to_contig_log.log
#        # parse the "realign.sam", and pop the perfect match to the 1st_round_contigs to generate the alternative remain_P1.fasta and remain_P2.fasta
#        python3 parse_contig_realign.py -fs ${haplotype_sam_path}/${allele_name}_realign_${cluster_id}.sam \
#                                        -fo ${haplotype_sam_path}/${allele_name}_remain_${cluster_id}.fasta \
#                                        -it_region ${flank_region} \
#                                        --contig_file ${contig_path}/${allele_name}_contig_${cluster_id}.fasta \
#                                        --allele_file ${raw_seq_path}/${allele_name}_${cluster_id}_allele.fasta \
#                                        --corrected_contig_output_file ${flanking_result_path}/corrected_contigs.fasta \
#                                        > ${haplotype_sam_path}/${allele_name}_realign_${cluster_id}.rpt
#    fi
#done
echo "[AIRRCall] [FLANKING SEQUENCE] Finish flanking sequence calling."

