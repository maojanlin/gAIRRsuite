# the out most directory
outer_dir=$1

# output reads and alleles fasta according to clusters
mkdir $outer_dir"asm_reads_and_alleles/"
#python3 output_reads.py -fa ../plot_tree/all_probs.fasta  -fr NA12878_S46_merged.fasta -fnp all_top600_allelelen_tcrv_idthrd100.cluster.pickle -fo $outer_dir"asm_reads_and_alleles/TCRV.fasta"
python3 output_reads.py -fa ../plot_tree/all_probs.fasta  -fr NA12878_S46_merged.fasta -fnp NA12878_tcrv_support.cluster.pickle -fo $outer_dir"asm_reads_and_alleles/TCRV.fasta"

# call the assembler spades.py
for cluster_id in {0..307}
do
    python3 ../SPAdes-3.11.1-Darwin/bin/spades.py -1 $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_H1.fasta" -2 $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_H2.fasta" --only-assembler -t 8 -o $outer_dir"asm_contigs/TCRV_"$cluster_id
    cp $outer_dir"asm_contigs/TCRV_"$cluster_id"/contigs.fasta" $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta"
done

# align alleles.fasta back to contigs.fasta
mkdir $outer_dir"asm_alignment/"
for cluster_id in {0..307}
do
    ../bwa/bwa index     $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta"
    ../bwa/bwa mem -t 16 $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta" $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_allele.fasta" > $outer_dir"asm_alignment/align_TCRV_"$cluster_id".sam"
done
