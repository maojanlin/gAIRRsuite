# output reads and alleles fasta according to clusters
python3 output_reads.py -fa ../plot_tree/all_probs.fasta  -fr NA12878_S46_merged.fasta -fnp all_top600_allelelen_ig_idthrd100.cluster.pickle -fo asm_reads_and_alleles/IG.fasta

# call the assembler spades.py
for cluster_id in {0..62}
do
    python3 ../SPAdes-3.11.1-Darwin/bin/spades.py -1 "./asm_reads_and_alleles/IG_"$cluster_id"_H1.fasta" -2 "./asm_reads_and_alleles/IG_"$cluster_id"_H2.fasta" --only-assembler -t 8 -o "./asm_contigs/IG_"$cluster_id
    cp "./asm_contigs/IG_"$cluster_id"/contigs.fasta" "./asm_contigs/IG_contigs_"$cluster_id".fasta"
done

# align alleles.fasta to contigs.fasta
for cluster_id in {0..62}
do
    ../bwa/bwa index "./asm_contigs/IG_contigs_"$cluster_id".fasta"
    ../bwa/bwa mem "./asm_contigs/IG_contigs_"$cluster_id".fasta" "./asm_reads_and_alleles/IG_"$cluster_id"_allele.fasta" > "./asm_alignment/align_IG_"$cluster_id".sam"
done

# record the results
touch ./asm_alignment/assembly_call.txt
for cluster_id in {0..62}
do
    python3 parse_bwa_sam.py -td 0 -fs "./asm_alignment/align_IG_"$cluster_id".sam" >> ./asm_alignment/assembly_call.txt
done

