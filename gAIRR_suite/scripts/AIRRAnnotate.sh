# the out most directory
outer_dir="target_annotation/"
list_allele_name="TCRV TCRJ BCRV BCRJ TCRD_plusHep BCRD_plusHep"
allele_dir="./gAIRR_suite/material/"
allele_suffix="_alleles_parsed.fasta"


person_name="HG002-part"
asm_path_H1="./gAIRR_suite/example/samples/HG002-S22-H1-000000F_1900000-2900000.fasta"

# setting for the data
echo "[AIRRAnnotate] Indexing assembly..."
bwa index ${asm_path_H1}
#bwa index ${asm_path_H2}

mkdir -p ${outer_dir}
for allele_name in ${list_allele_name}; do
    echo "[AIRRAnnotate] Align IMGT ${allele_name} alleles to assembly..."
    if [ ${allele_name} == "TCRD_plusHep" ] || [ ${allele_name} == "BCRD_plusHep" ]; then
        echo "[AIRRAnnotate] Adjust BWA parameters for shorter alleles..."
        bwa mem -t 16 -a -T 10 ${asm_path_H1} ${allele_dir}${allele_name}${allele_suffix} > ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H1.sam
#        bwa mem -t 16 -a -T 10 ${asm_path_H2} ${allele_dir}${allele_name}${allele_suffix} > ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H2.sam
    else
        bwa mem -t 16 -a ${asm_path_H1} ${allele_dir}${allele_name}${allele_suffix} > ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H1.sam
#        bwa mem -t 16 -a ${asm_path_H2} ${allele_dir}${allele_name}${allele_suffix} > ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H2.sam
    fi


    echo "[AIRRAnnotate] Parse the ${allele_name} alleles sam files to annotation report..."
    python3 gAIRR_suite/scripts/annotation_with_asm.py -foa   ${outer_dir}annotation_${person_name}_${allele_name}.txt \
                                                       -foma  ${outer_dir}annotation_imperfect_${person_name}_${allele_name}.txt \
                                                       -fom   ${outer_dir}novel_${person_name}_${allele_name}.fasta \
                                                       -fof   ${outer_dir}flanking_${person_name}_${allele_name}.fasta \
                                                       -fos   ${outer_dir}summary_${person_name}_${allele_name}.rpt \
                                                       -ext   200 \
                                                       -fs1   ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H1.sam \
                                                       -fp1   ${outer_dir}dict_${person_name}_asm_H1.pickle \
                                                       -fasm1 ${asm_path_H1} \
#                                                       -fs2   ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H2.sam \
#                                                       -fp2   ${outer_dir}dict_${person_name}_asm_H2.pickle \
#                                                       -fasm2 ${asm_path_H2}
done
echo "[ANNOTATION WITH ASM] Finished!"
