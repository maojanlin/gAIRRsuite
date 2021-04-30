# the out most directory
outer_dir="RSS_checking/"
report_dir="report/"

hep_nona_file="TR_heptamer_nonamer.fasta"
person_name="PL"
person_name="HG007"
list_allele="TCRV TCRJ"

# setting for the data
echo "[AIRRCall] Checking RSS in Flanking Regions..."

mkdir -p ${outer_dir}
mkdir -p ${outer_dir}/${report_dir}
for allele in ${list_allele}; do
    flanking_file="RSS_checking/${person_name}_${allele}_flanking_haplotypes.fasta"
    bwa index ${flanking_file}
    
    bwa mem -t 16 -a -T 7 -k 7 ${flanking_file} ${hep_nona_file} > ${outer_dir}bwa_${person_name}_${allele}_RSS_to_flanking.sam
    samtools sort ${outer_dir}bwa_${person_name}_${allele}_RSS_to_flanking.sam > ${outer_dir}bwa_${person_name}_${allele}_RSS_to_flanking_sorted.bam
    samtools index ${outer_dir}bwa_${person_name}_${allele}_RSS_to_flanking_sorted.bam
    
    python3 check_RSS.py -fs  ${outer_dir}/bwa_${person_name}_${allele}_RSS_to_flanking.sam \
                         -gtype ${allele} \
                         -fo  ${outer_dir}/${report_dir}/${person_name}_${allele}_RSS.csv \
                         -fos ${outer_dir}/${report_dir}/summary_${person_name}_${allele}_RSS.rpt
done
echo "[AIRRCall] RSS Checking Finished!"
