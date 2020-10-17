# the out most directory
outer_dir="target_annotation/"
list_allele_name="TCRV TCRJ BCRV BCRJ TCRD_plusHep BCRD_plusHep"
allele_dir="../genomeData/"
allele_suffix="_alleles_parsed.fasta"
#person_name="NA24385"
#asm_path_H1="../../asm/NA24385/HG002-H1.fa"
#asm_path_H2="../../asm/NA24385/HG002-H2.fa"
#person_name="CHM13-T2T"
#asm_path_H1="../genomeData/T2T_CHM13/chm13.draft_v1.0.fasta"
person_name="GCA_HG004"
asm_path_H1="../genomeData/pacbio_GIAB/GCA_001549595.1_GIAB_Ashkenazim_Mother_HG004_NA24143_hu8E87A9_PacBio_Assembly_with_PBcR_genomic.fna"
#asm_path_H1="../genomeData/pacbio_GIAB/GCA_001549605.1_GIAB_Ashkenazim_Father_HG003_NA24149_hu6E4515_PacBio_Assembly_with_PBcR_genomic.fna"

# setting for the data
echo "[AIRRAnnotate] Indexing assembly..."
#bwa index ${asm_path_H1}
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
    python3 scripts/annotation_with_asm.py -foa   ${outer_dir}annotation_${person_name}_${allele_name}.txt \
                                           -foma  ${outer_dir}annotation_imperfect_${person_name}_${allele_name}.txt \
                                           -fom   ${outer_dir}novel_${person_name}_${allele_name}.fasta \
                                           -fof   ${outer_dir}flanking_${person_name}_${allele_name}.fasta \
                                           -fos   ${outer_dir}summary_${person_name}_${allele_name}.rpt \
                                           -ext   200 \
                                           -fs1   ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H1.sam \
                                           -fp1   ${outer_dir}dict_${person_name}_asm_H1.pickle \
                                           -fasm1 ${asm_path_H1} \
#                                           -fs2   ${outer_dir}bwa_${person_name}_${allele_name}_alleles_to_asm_H2.sam \
#                                           -fp2   ${outer_dir}dict_${person_name}_asm_H2.pickle \
#                                           -fasm2 ${asm_path_H2}
done
echo "[ANNOTATION WITH ASM] Finished!"
