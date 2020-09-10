# path parameters
outer_dir=${1}/${4}_${2}_novel/ 
allele_name=$2
allele_path=$3  
read_path_1=$5  
read_path_2=$6  

# setting for the data
mkdir -p ${outer_dir}
#bwa index ${allele_path}
#bwa mem -t 16 ${allele_path} ${read_path_1} ${read_path_2} > ${outer_dir}bwa_read_to_allele.sam

# start analysis
echo "[AIRRCall] [NOVEL ALLELE] Finding novel alleles..."
rm ${outer_dir}corrected_alleles_raw.fasta
python3 parse_cluster_realign.py -fs  ${outer_dir}bwa_read_to_allele.sam \
                                 -fc  ${allele_path} \
                                 -for test.txt \
                                 -foc ${outer_dir}corrected_alleles_raw.fasta \
                                 >    ${outer_dir}bwa_alleles_cover_link.rpt \
                                 2>   ${outer_dir}bwa_alleles_all.txt

# product refinement, bwa mem -a mode is used to find all duplications
echo "[AIRRCall] [NOVEL ALLELE] Product refinement..."
python3 filter_corrected_alleles.py -fa  ${allele_path} \
                                    -fca ${outer_dir}corrected_alleles_raw.fasta \
                                    -fof ${outer_dir}corrected_alleles_filtered.fasta \
                                    -foe ${outer_dir}corrected_alleles_extended.fasta

# merge the refined product with original alleles
echo "[AIRRCall] [NOVEL ALLELE] Add novel allele to target alleles..."
python3 merge_novel_alleles.py -fa  ${allele_path} \
                               -fn  ${outer_dir}corrected_alleles_filtered.fasta \
                               -fom ${outer_dir}/${allele_name}_with_novel.fasta

echo "[AIRRCall] [NOVEL ALLELE] Finished!"
