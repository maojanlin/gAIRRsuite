# path of AIRRAnnotate
asm_path_H1="../../asm/NA24385/HG002-H1.fa"
asm_path_H2="../../asm/NA24385/HG002-H2.fa"
annotation_dir="target_annotation/"

# path of AIRRCall
workspace="target_call"
allele_name="BCRV"
person_name="NA24385"


# verification
novel_dir=${workspace}/${person_name}_${allele_name}_novel/
echo "[AIRRVerify] Verify novel allele calls..."
mkdir -p ${novel_dir}verification/
bwa mem -t 16 ${asm_path_H1} ${novel_dir}corrected_alleles_filtered.fasta > ${novel_dir}verification/bwa_c_alleles_f_verify_H1.sam 
bwa mem -t 16 ${asm_path_H2} ${novel_dir}corrected_alleles_filtered.fasta > ${novel_dir}verification/bwa_c_alleles_f_verify_H2.sam 
python3 scripts/verify_novel_alleles.py -fs1 ${novel_dir}verification/bwa_c_alleles_f_verify_H1.sam \
                                        -fs2 ${novel_dir}verification/bwa_c_alleles_f_verify_H2.sam \
                                        -fo  ${novel_dir}verification/novel_verification.rpt


flank_dir=${workspace}/${person_name}_${allele_name}_flanking/
echo "[AIRRVerify] Verify flanking sequences..."
bwa mem -t 16 -a ${asm_path_H1} ${flank_dir}/flanking_result/flanking_haplotypes.fasta > ${flank_dir}/flanking_result/bwa_flanking_haplotype_H1.sam
bwa mem -t 16 -a ${asm_path_H2} ${flank_dir}/flanking_result/flanking_haplotypes.fasta > ${flank_dir}/flanking_result/bwa_flanking_haplotype_H2.sam

mkdir -p ${flank_dir}verification/
python3 scripts/annotation_with_asm.py -foa  ${flank_dir}verification/annotation_${allele_name}.txt \
                                       -foma ${flank_dir}verification/annotation_imperfect_${allele_name}.txt \
                                       -fom  ${flank_dir}verification/corrected_${allele_name}_flanking.fasta \
                                       -fs1  ${flank_dir}flanking_result/bwa_flanking_haplotype_H1.sam \
                                       -fs2  ${flank_dir}flanking_result/bwa_flanking_haplotype_H2.sam

python3 scripts/compare_annotation.py  -fab ${annotation_dir}annotation_imperfect_${person_name}_${allele_name}.txt \
                                       -fan ${flank_dir}verification/annotation_imperfect_${allele_name}.txt \
                                       -fh  ${flank_dir}flanking_result/flanking_haplotypes.fasta \
                                       -for ${flank_dir}verification/annotation_report.txt \
                                       -fos ${flank_dir}verification/annotation_summary.txt

