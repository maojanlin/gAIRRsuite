# the out most directory
outer_dir="merged_IMGT_alleles/"
allele_path="../plot_tree/BCRV_alleles.fasta"
novel_path="CHM13_bcrv_novel_alleles/corrected_alleles_filtered.fasta"
output_name="BCRV_with_CHM13_novel_alleles.fasta"


echo "[MERGE NOVEL] Merge novel alleles"
mkdir -p ${outer_dir}
python3 merge_novel_alleles.py -fa  ${allele_path} \
                               -fn  ${novel_path} \
                               -fom ${outer_dir}${output_name}
echo "[MERGE NOVEL] Finished!"
