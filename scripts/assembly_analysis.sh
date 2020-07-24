# the out most directory
outer_dir="NA12878_tcrv_support_asm/"

# building the environment and do the 1st round assembly
#./assembly_1st_round.sh $outer_dir

# record the results
rm $outer_dir"asm_alignment/assembly_call.txt"
rm $outer_dir"asm_alignment/flanking_contigs.fasta"
rm $outer_dir"asm_alignment/flanking_size300.fasta"
#for cluster_id in 106
for cluster_id in {0..307}
do
    echo $cluster_id
    flank_region=$(python3 parse_bwa_sam.py -fs $outer_dir"asm_alignment/align_TCRV_"$cluster_id".sam" -fc $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta" -td 0 -fo $outer_dir"asm_alignment/assembly_call.txt" -fr $outer_dir"asm_alignment/flanking_contigs.fasta" --cluster_id $cluster_id"-1" -fsize 300 -frs $outer_dir"asm_alignment/flanking_size300.fasta")
    # if parse_bwa_sam.py report region, then do the 2nd round assembly
    if [ ! -z $flank_region ]  # if $flank_region exist
    then
        ./assembly_2nd_round.sh $outer_dir $cluster_id $flank_region $flank_region
        python3 parse_bwa_sam.py -fs $outer_dir"asm_2nd_alignment/align_TCRV_2nd_"$cluster_id".sam" -fc $outer_dir"asm_2nd_contigs/TCRV_2nd_contigs_"$cluster_id".fasta" -td 0 -fo $outer_dir"asm_alignment/assembly_call.txt" -fr $outer_dir"asm_alignment/flanking_contigs.fasta" --cluster_id $cluster_id"-2" -fsize 300 -frs $outer_dir"asm_alignment/flanking_size300.fasta"
    fi
done

