
asm_path_H1="../asm_NA12878/NA12878-H1.fa"
asm_path_H2="../asm_NA12878/NA12878-H2.fa"
outer_dir="target_call/NA12878_TCRJ_novel/"


# verification
echo "[AIRRVerify] Verify novel allele calls..."
mkdir -p ${outer_dir}verification/
../bwa/bwa mem -t 16 ${asm_path_H1} ${outer_dir}corrected_alleles_filtered.fasta > ${outer_dir}verification/bwa_c_alleles_f_verify_H1.sam 
../bwa/bwa mem -t 16 ${asm_path_H2} ${outer_dir}corrected_alleles_filtered.fasta > ${outer_dir}verification/bwa_c_alleles_f_verify_H2.sam 
python3 verify_novel_alleles.py -fs1 ${outer_dir}verification/bwa_c_alleles_f_verify_H1.sam \
                                -fs2 ${outer_dir}verification/bwa_c_alleles_f_verify_H2.sam \
                                -fo  ${outer_dir}verification/verification.rpt

