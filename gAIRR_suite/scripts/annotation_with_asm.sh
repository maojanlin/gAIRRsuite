# the out most directory
outer_dir="asm_annotation/"
allele_path_tcrv="../plot_tree/TCRV_alleles.fasta"
allele_path_tcrj="../plot_tree/TCRJ_alleles.fasta"
allele_path_bcrv="../plot_tree/BCRV_alleles.fasta"
asm_path_H1="../asm_NA12878/NA12878-H1.fa"
asm_path_H2="../asm_NA12878/NA12878-H2.fa"

# setting for the data
echo "[ANNOTATION WITH ASM] Align IMGT alleles to assembly..."
mkdir -p ${outer_dir}
../bwa/bwa index ${asm_path_H1}
../bwa/bwa index ${asm_path_H2}
../bwa/bwa mem -t 16 -a ${asm_path_H1} ${allele_path_tcrv} > ${outer_dir}bwa_NA12878_tcrv_alleles_to_asm_H1.sam
../bwa/bwa mem -t 16 -a ${asm_path_H2} ${allele_path_tcrv} > ${outer_dir}bwa_NA12878_tcrv_alleles_to_asm_H2.sam
../bwa/bwa mem -t 16 -a ${asm_path_H1} ${allele_path_tcrj} > ${outer_dir}bwa_NA12878_tcrj_alleles_to_asm_H1.sam
../bwa/bwa mem -t 16 -a ${asm_path_H2} ${allele_path_tcrj} > ${outer_dir}bwa_NA12878_tcrj_alleles_to_asm_H2.sam
../bwa/bwa mem -t 16 -a ${asm_path_H1} ${allele_path_bcrv} > ${outer_dir}bwa_NA12878_bcrv_alleles_to_asm_H1.sam
../bwa/bwa mem -t 16 -a ${asm_path_H2} ${allele_path_bcrv} > ${outer_dir}bwa_NA12878_bcrv_alleles_to_asm_H2.sam

# start analysis
echo "[ANNOTATION WITH ASM] Parse the sam files to annotation report..."
echo "[TCRV Alleles]"
python3 annotation_with_asm.py -foa ${outer_dir}annotation_NA12878_tcrv.txt \
                               -foma ${outer_dir}annotation_imperfect_NA12878_tcrv.txt \
                               -fom ${outer_dir}corrected_NA12878_tcrv.fasta \
                               -fs1 ${outer_dir}bwa_NA12878_tcrv_alleles_to_asm_H1.sam \
                               -fs2 ${outer_dir}bwa_NA12878_tcrv_alleles_to_asm_H2.sam
echo "[TCRJ Alleles]"
python3 annotation_with_asm.py -foa ${outer_dir}annotation_NA12878_tcrj.txt \
                               -foma ${outer_dir}annotation_imperfect_NA12878_tcrj.txt \
                               -fom ${outer_dir}corrected_NA12878_tcrj.fasta \
                               -fs1 ${outer_dir}bwa_NA12878_tcrj_alleles_to_asm_H1.sam \
                               -fs2 ${outer_dir}bwa_NA12878_tcrj_alleles_to_asm_H2.sam
echo "[BCRV Alleles]"
python3 annotation_with_asm.py -foa ${outer_dir}annotation_NA12878_bcrv.txt \
                               -foma ${outer_dir}annotation_imperfect_NA12878_bcrv.txt \
                               -fom ${outer_dir}corrected_NA12878_bcrv.fasta \
                               -fs1 ${outer_dir}bwa_NA12878_bcrv_alleles_to_asm_H1.sam \
                               -fs2 ${outer_dir}bwa_NA12878_bcrv_alleles_to_asm_H2.sam
echo "[ANNOTATION WITH ASM] Finished!"
