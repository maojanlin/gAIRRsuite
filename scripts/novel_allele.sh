# the out most directory
outer_dir="NA12878_tcrv_novel_alleles/"
allele_path="../plot_tree/TCRV_alleles.fasta"
read_path_1="NA12878_S46_R1.fasta"
read_path_2="NA12878_S46_R2.fasta"
asm_path_H1="../asm_NA12878/NA12878-H1.fa"
asm_path_H2="../asm_NA12878/NA12878-H2.fa"

# setting for the data
mkdir -p ${outer_dir}
../bwa/bwa index ${allele_path}
../bwa/bwa mem -t 16 ${allele_path} ${read_path_1} ${read_path_2} > ${outer_dir}bwa_read_to_allele.sam

# start analysis
echo "[NOVEL ALLELE] Finding novel alleles..."
rm ${outer_dir}corrected_alleles_raw.fasta
python3 parse_cluster_realign.py -fs  ${outer_dir}bwa_read_to_allele.sam \
                                 -fc  ${allele_path} \
                                 -for test.txt \
                                 -foc ${outer_dir}corrected_alleles_raw.fasta \
                                 >    ${outer_dir}bwa_alleles_cover_link.rpt \
                                 2>   ${outer_dir}bwa_alleles_all.txt

# product refinement, bwa mem -a mode is used to find all duplications
echo "[NOVEL ALLELE] Product refinement..."
../bwa/bwa index ${outer_dir}corrected_alleles_raw.fasta
../bwa/bwa mem -t 16 -a ${outer_dir}corrected_alleles_raw.fasta ${allele_path} > ${outer_dir}bwa_c_alleles_realign.sam
python3 filter_corrected_alleles.py -fs  ${outer_dir}bwa_c_alleles_realign.sam \
                                    -fca ${outer_dir}corrected_alleles_raw.fasta \
                                    -fof ${outer_dir}corrected_alleles_filtered.fasta

# verification
echo "[NOVEL ALLELE] Verify the product..."
mkdir -p ${outer_dir}verification/
../bwa/bwa mem -t 16 ${asm_path_H1} ${outer_dir}corrected_alleles_filtered.fasta > ${outer_dir}verification/bwa_c_alleles_f_verify_H1.sam 
../bwa/bwa mem -t 16 ${asm_path_H2} ${outer_dir}corrected_alleles_filtered.fasta > ${outer_dir}verification/bwa_c_alleles_f_verify_H2.sam 
python3 verify_novel_alleles.py -fs1 ${outer_dir}verification/bwa_c_alleles_f_verify_H1.sam \
                                -fs2 ${outer_dir}verification/bwa_c_alleles_f_verify_H2.sam \
                                -fo  ${outer_dir}verification/verification.rpt

echo "[NOVEL ALLELE] Finished!"
