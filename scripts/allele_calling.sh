#annotation_path="asm_annotation/annotation_NA12878_tcrv.txt"
# path parameters
outer_dir=$1/$4_$2/
allele_path=$3
read_path_1=$5
read_path_2=$6
coverage_thrsd=100

# setting for the data
mkdir -p ${outer_dir}
bwa index ${allele_path}
bwa mem -t 16 -a ${allele_path} ${read_path_1} ${read_path_2} > ${outer_dir}bwa_read_to_allele_all.sam

# start analysis
echo "[AIRRCall] [ALLELE CALLING] Calling alleles..."
python3 scripts/analyze_read_depth_with_bwa.py -fs  ${outer_dir}bwa_read_to_allele_all.sam \
                                               -fa  ${allele_path} \
                                               -t   ${coverage_thrsd} \
                                               -foc ${outer_dir}read_depth_calling_by_bwa.rpt \
                                               -fop ${outer_dir}allele_support_reads.pickle \
                                               #-fv  ${annotation_path}
echo "[AIRRCall] [ALLELE CALLING] Finished!"
