# the out most directory
outer_dir="NA24385_tcrv_allele_with_novel_calling/"
#allele_path="../plot_tree/TCRJ_alleles.fasta"
allele_path="merged_IMGT_alleles/TCRV_with_NA24385_novel_alleles.fasta"
read_path_1="NA24385_S42_R1.fasta"
read_path_2="NA24385_S42_R2.fasta"
annotation_path="asm_annotation/annotation_NA24385_tcrv.txt"
coverage_thrsd=100

# setting for the data
mkdir -p ${outer_dir}
../bwa/bwa index ${allele_path}
../bwa/bwa mem -t 16 -a ${allele_path} ${read_path_1} ${read_path_2} > ${outer_dir}bwa_read_to_allele_all.sam

# start analysis
echo "[ALLELE CALLING] Calling alleles..."
python3 analyze_read_depth_with_bwa.py -fs  ${outer_dir}bwa_read_to_allele_all.sam \
                                       -fa  ${allele_path} \
                                       -t   ${coverage_thrsd} \
                                       -foc ${outer_dir}read_depth_calling_by_bwa.rpt \
                                       -fv  ${annotation_path}
echo "[ALLELE CALLING] Finished!"
