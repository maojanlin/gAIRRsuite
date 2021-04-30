# the out most directory
outer_dir="RSS_checking/"
first_dir="first_scan/"
second_dir="second_scan/"
report_dir="report/"

RSS_file="TR_RSS.fasta"
hep_nona_file="TR_heptamer_nonamer.fasta"
person_name="HQ-12"
person_name="HGSVC-24"
#person_name="PL"
list_allele="TCRV TCRJ"

# setting for the data
echo "[AIRRCall] Checking RSS in Flanking Regions..."

mkdir -p ${outer_dir}
mkdir -p ${outer_dir}/${first_dir}
mkdir -p ${outer_dir}/${second_dir}
mkdir -p ${outer_dir}/${report_dir}
for allele in ${list_allele}; do
    flanking_file="RSS_checking/${person_name}_${allele}_flanking_haplotypes.fasta"
    bwa index ${flanking_file}
    
    bwa mem -t 16 -a -k 10 -T 20 ${flanking_file} ${RSS_file} > ${outer_dir}${first_dir}/bwa_${person_name}_${allele}_RSS_first_scan.sam
    #samtools sort  ${outer_dir}${first_dir}/bwa_${person_name}_${allele}_RSS_first_scan.sam > ${outer_dir}${first_dir}bwa_${person_name}_${allele}_RSS_first_scan_sorted.bam
    #samtools index ${outer_dir}${first_dir}/bwa_${person_name}_${allele}_RSS_first_scan_sorted.bam

    python3 check_RSS_1st_scan.py -fs   ${outer_dir}${first_dir}/bwa_${person_name}_${allele}_RSS_first_scan.sam \
                                  -ff   ${flanking_file} \
                                  -fod  ${outer_dir}${first_dir}/detailed_${person_name}_${allele}_first_scan.rpt \
                                  -fos  ${outer_dir}${first_dir}/summary_${person_name}_${allele}_first_scan.rpt \
                                  -fonf ${outer_dir}${first_dir}/novel_RSS_${person_name}_${allele}_first_scan.fasta \
                                  -fomf ${outer_dir}${first_dir}/missing_RSS_${person_name}_${allele}_first_scan.fasta

    bwa index ${outer_dir}${first_dir}/missing_RSS_${person_name}_${allele}_first_scan.fasta
    bwa mem -t 16 -a -T 7 -k 7 ${outer_dir}${first_dir}/missing_RSS_${person_name}_${allele}_first_scan.fasta ${hep_nona_file} > \
                               ${outer_dir}${second_dir}/bwa_${person_name}_${allele}_RSS_to_missing.sam
    #samtools sort ${outer_dir}${second_dir}/bwa_${person_name}_${allele}_RSS_to_missing.sam > ${outer_dir}${second_dir}bwa_${person_name}_${allele}_RSS_to_missing_sorted.bam
    #samtools index ${outer_dir}${second_dir}/bwa_${person_name}_${allele}_RSS_to_missing_sorted.bam
    
    python3 check_RSS_2nd_scan.py -fs    ${outer_dir}${second_dir}/bwa_${person_name}_${allele}_RSS_to_missing.sam \
                                  -ff    ${outer_dir}${first_dir}/missing_RSS_${person_name}_${allele}_first_scan.fasta \
                                  -gtype ${allele} \
                                  -fo    ${outer_dir}/${second_dir}/${person_name}_${allele}_RSS.csv \
                                  -fos   ${outer_dir}/${second_dir}/summary_${person_name}_${allele}_RSS.rpt \
                                  -fon   ${outer_dir}/${second_dir}/novel_RSS_${person_name}_${allele}_second_scan.fasta
done
echo "[AIRRCall] RSS First Scanning Finished!"
