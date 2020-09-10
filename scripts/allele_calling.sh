#annotation_path="asm_annotation/annotation_NA12878_tcrv.txt"
# path parameters
outer_dir=$1/$4_$2/   #"NA12878_tcrv_novel_alleles/"
allele_path=$3  #"../plot_tree/TCRV_alleles.fasta"
read_path_1=$5  #"NA12878_S46_R1.fasta"
read_path_2=$6  #"NA12878_S46_R2.fasta"
coverage_thrsd=100

# setting for the data
mkdir -p ${outer_dir}
bwa index ${allele_path}
bwa mem -t 16 -a ${allele_path} ${read_path_1} ${read_path_2} > ${outer_dir}bwa_read_to_allele_all.sam

# start analysis
echo "[AIRRCall] [ALLELE CALLING] Calling alleles..."
python3 analyze_read_depth_with_bwa.py -fs  ${outer_dir}bwa_read_to_allele_all.sam \
                                       -fa  ${allele_path} \
                                       -t   ${coverage_thrsd} \
                                       -foc ${outer_dir}read_depth_calling_by_bwa.rpt \
                                       -fop ${outer_dir}allele_support_reads.pickle \
                                       #-fv  ${annotation_path}
echo "[AIRRCall] [ALLELE CALLING] Finished!"
