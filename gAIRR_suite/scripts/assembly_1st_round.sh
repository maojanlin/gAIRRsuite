# the out most directory
outer_dir=$1
human_id=$2

# output reads and alleles fasta according to clusters
mkdir $outer_dir"asm_reads_and_alleles/"
#python3 output_reads.py -fa ../plot_tree/all_probs.fasta  -fr NA12878_S46_merged.fasta -fnp all_top600_allelelen_tcrv_idthrd100.cluster.pickle -fo $outer_dir"asm_reads_and_alleles/TCRV.fasta"
#python3 output_reads.py -fa ../plot_tree/all_probs.fasta  -fr $human_id"_S42_merged.fasta" -fnp $human_id"_tcrv_support.cluster.pickle" -fo $outer_dir"asm_reads_and_alleles/TCRV.fasta"
python3 ../../immunogenomics/scripts/output_reads.py -fr $human_id"_S42_merged.fasta" -fa clustering/TR_all.fasta -fnp $human_id"/"$human_id"_tcrv_support.cluster.pickle" -fo $outer_dir"asm_reads_and_alleles/TCRV.fasta"

# call the assembler spades.py
for cluster_id in {0..307}
do
    python3 ../../SPAdes-3.11.1-Linux/bin/spades.py -1 $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_H1.fasta" -2 $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_H2.fasta" --only-assembler -t 8 -o $outer_dir"asm_contigs/TCRV_"$cluster_id
    python3 ../../immunogenomics/scripts/copy_1st_contig.py -fc $outer_dir"asm_contigs/TCRV_"$cluster_id"/contigs.fasta" -fo $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta"
done

# align alleles.fasta back to contigs.fasta
mkdir $outer_dir"asm_alignment/"
for cluster_id in {0..307}
do
    bwa index     $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta"
    bwa mem -t 16 $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta" $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_allele.fasta" > $outer_dir"asm_alignment/align_TCRV_"$cluster_id".sam"
done
