# the out most directory
outer_dir="target_annotation/"
list_allele_name="TCRV TCRJ BCRV"
allele_dir="../genomeData/"
allele_suffix="_alleles_parsed.fasta"
person_name="NA24385"
asm_path_H1="../../asm/NA24385/HG002-H1.fa"
asm_path_H2="../../asm/NA24385/HG002-H2.fa"

# setting for the data
echo "[AIRRAnnotate] Indexing assembly..."
#bwa index ${asm_path_H1}
#bwa index ${asm_path_H2}

mkdir -p ${outer_dir}
for allele_name in ${list_allele_name}; do
    echo "[AIRRAnnotate] Align IMGT ${allele_name} alleles to assembly..."
    bwa mem -t 16 -a ${asm_path_H1} ${allele_dir}${allele_name}${allele_suffix} > ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H1.sam
    bwa mem -t 16 -a ${asm_path_H2} ${allele_dir}${allele_name}${allele_suffix} > ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H2.sam

    echo "[AIRRAnnotate] Parse the ${allele_name} alleles sam files to annotation report..."
    python3 scripts/annotation_with_asm.py -foa   ${outer_dir}annotation_${person_name}_${allele_name}.txt \
                                           -foma  ${outer_dir}annotation_imperfect_${person_name}_${allele_name}.txt \
                                           -fom   ${outer_dir}novel_${person_name}_${allele_name}.fasta \
                                           -fof   ${outer_dir}flanking_${person_name}_${allele_name}.fasta \
                                           -ext   200 \
                                           -fs1   ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H1.sam \
                                           -fasm1 ${asm_path_H1} \
                                           -fs2   ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H2.sam \
                                           -fasm2 ${asm_path_H2}
done
echo "[ANNOTATION WITH ASM] Finished!"
