#annotation_path="asm_annotation/annotation_NA12878_tcrv.txt"
# path parameters
outer_dir=$1/$4_$2/
allele_name=$2
allele_path=$3
read_path_1=$5
read_path_2=$6
thread=$7
coverage_thrsd=100

# setting for the data
mkdir -p ${outer_dir}
bwa index ${allele_path}
if [ ${allele_name} == "TCRD_plusHep" ] || [ ${allele_name} == "BCRD_plusHep" ]; then
    echo "[AIRRCall] Adjust BWA parameters for shorter alleles..."
    bwa mem -t ${thread} -a -T 10 ${allele_path} ${read_path_1} ${read_path_2} > ${outer_dir}bwa_read_to_allele_all.sam
else
    bwa mem -t ${thread} -a ${allele_path} ${read_path_1} ${read_path_2} > ${outer_dir}bwa_read_to_allele_all.sam
fi

# start analysis
echo "[AIRRCall] [ALLELE CALLING] Calling alleles..."
python3 scripts/analyze_read_depth_with_bwa.py -fs  ${outer_dir}bwa_read_to_allele_all.sam \
                                               -fa  ${allele_path} \
                                               -t   ${coverage_thrsd} \
                                               -foc ${outer_dir}read_depth_calling_by_bwa.rpt \
                                               -fop ${outer_dir}allele_support_reads.pickle \
                                               #-fv  ${annotation_path}

python3 scripts/calling_threshold.py -dp ${outer_dir}read_depth_calling_by_bwa.rpt > ${outer_dir}gAIRR-call_report.rpt
echo "[AIRRCall] [ALLELE CALLING] Finished!"
