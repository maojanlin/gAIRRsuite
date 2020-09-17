# the out most directory
outer_dir="target_annotation/"
list_allele_name="TCRV TCRJ BCRV BCRJ TCRD_plusHep BCRD_plusHep"
novel_prefix="novel_"
novel_suffix=".fasta"
flank_prefix="flanking_"
flank_suffix=".fasta"

mkdir -p ${outer_dir}
for allele_name in ${list_allele_name}; do
    echo "[Database Collect] Merge ${allele_name} novel alleles..."
    python3 ./scripts/collect_novel_alleles.py -fa "$( ls ${outer_dir}/${novel_prefix}*${allele_name}${novel_suffix} )" \
                                               -for ${outer_dir}/database_novel_${allele_name}.rpt \
                                               -fof ${outer_dir}/database_novel_${allele_name}.fasta
    
    echo "[Database Collect] Merge ${allele_name} flanking sequences..."
    python3 ./scripts/collect_novel_alleles.py -fa "$( ls ${outer_dir}/${flank_prefix}*${allele_name}${novel_suffix} )" \
                                               -for ${outer_dir}/database_flanking_${allele_name}.rpt \
                                               -fof ${outer_dir}/database_flanking_${allele_name}.fasta
done
echo "[Database Collect] Finished!"
