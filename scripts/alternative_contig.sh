# the out most directory
outer_dir="NA12878_tcrv_support_asm/"
cluster_id="106"

../bwa/bwa mem $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta" $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_H1.fasta" $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_H2.fasta" > $outer_dir"asm_realign/TCRV_realign_"$cluster_id".sam"
python3 parse_contig_realign.py -fs $outer_dir"asm_realign/TCRV_realign_"$cluster_id".sam" -fo $outer_dir"asm_realign/TCRV_remain_"$cluster_id".fasta" > "TCRV_realign_"$cluster_id".rpt"
# align remain.fasta back to contig
../bwa/bwa mem $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta" $outer_dir"asm_realign/TCRV_remain_"$cluster_id"_P1.fasta" $outer_dir"asm_realign/TCRV_remain_"$cluster_id"_P2.fasta" > $outer_dir"asm_realign/align_TCRV_remain_"$cluster_id".sam"
#samtools sort $outer_dir"asm_realign/align_TCRV_remain_"$cluster_id".sam" > $outer_dir"asm_realign/align_TCRV_extreme_remain_"$cluster_id"_sorted.bam"
#samtools index $outer_dir"asm_realign/align_TCRV_remain_"$cluster_id"_sorted.bam"
python3 ../SPAdes-3.11.1-Darwin/bin/spades.py -1 $outer_dir"asm_realign/TCRV_remain_"$cluster_id"_P1.fasta" -2 $outer_dir"asm_realign/TCRV_remain_"$cluster_id"_P2.fasta" --only-assembler -t 8 -o $outer_dir"asm_realign/TCRV_remain_"$cluster_id


# align to assembly for verification
cp $outer_dir"asm_realign/TCRV_remain_"$cluster_id"/contigs.fasta" $outer_dir"asm_realign/TCRV_contigs_remain_"$cluster_id".fasta"
cat $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta" $outer_dir"asm_realign/TCRV_contigs_remain_"$cluster_id".fasta" > $outer_dir"asm_realign/TCRV_contigs_pair_"$cluster_id".fasta"
../bwa/bwa mem ../asm_NA12878/NA12878-H1.fa $outer_dir"asm_realign/TCRV_contigs_pair_"$cluster_id".fasta" > $outer_dir"asm_realign/align_pair_"$cluster_id"_H1.sam"
#samtools sort $outer_dir"asm_realign/align_pair_"$cluster_id"_H1.sam" > $outer_dir"asm_realign/align_pair_"$cluster_id"_H1_sorted.bam"
#samtools index $outer_dir"asm_realign/align_pair_"$cluster_id"_H1_sorted.bam"
../bwa/bwa mem ../asm_NA12878/NA12878-H2.fa $outer_dir"asm_realign/TCRV_contigs_pair_"$cluster_id".fasta" > $outer_dir"asm_realign/align_pair_"$cluster_id"_H2.sam"
#samtools sort $outer_dir"asm_realign/align_pair_"$cluster_id"_H2.sam" > $outer_dir"asm_realign/align_pair_"$cluster_id"_H2_sorted.bam"
#samtools index $outer_dir"asm_realign/align_pair_"$cluster_id"_H2_sorted.bam"


