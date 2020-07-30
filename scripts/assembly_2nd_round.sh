# the out most directory
outer_dir=$1
cluster_id=$2
interested_region=$3

# align reads back to the SPAdes_assembly contigs as "realign.sam"
../bwa/bwa mem -t 16 $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta" $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_H1.fasta" $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_H2.fasta" > $outer_dir"asm_back_align/TCRV_realign_"$cluster_id".sam"
# parse the "realign.sam", and pop the perfect match to the 1st_round_contigs to generate the alternative remain_P1.fasta and remain_P2.fasta
python3 parse_contig_realign.py -fs $outer_dir"asm_back_align/TCRV_realign_"$cluster_id".sam" -fo $outer_dir"asm_back_align/TCRV_remain_"$cluster_id".fasta" > $outer_dir"asm_back_align/TCRV_realign_"$cluster_id".rpt" -it_region $interested_region --contig_file $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta" --allele_file $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_allele.fasta" --corrected_contig_output_file $outer_dir"asm_alignment/corrected_contig_size300.fasta"

## align remain.fasta back to 1st_round_contig for evaluation
#../bwa/bwa mem -t 16 $outer_dir"asm_contigs/TCRV_contigs_"$cluster_id".fasta" $outer_dir"asm_back_align/TCRV_remain_"$cluster_id"_P1.fasta" $outer_dir"asm_back_align/TCRV_remain_"$cluster_id"_P2.fasta" > $outer_dir"asm_back_align/align_TCRV_remain_"$cluster_id".sam"
#samtools sort $outer_dir"asm_back_align/align_TCRV_remain_"$cluster_id".sam" > $outer_dir"asm_back_align/align_TCRV_extreme_remain_"$cluster_id"_sorted.bam"
#samtools index $outer_dir"asm_back_align/align_TCRV_remain_"$cluster_id"_sorted.bam"



# use SPAdes to assembly 2nd_round_contig
python3 ../SPAdes-3.11.1-Darwin/bin/spades.py -1 $outer_dir"asm_back_align/TCRV_remain_"$cluster_id"_P1.fasta" -2 $outer_dir"asm_back_align/TCRV_remain_"$cluster_id"_P2.fasta" --only-assembler -t 8 -o $outer_dir"asm_2nd_contigs/TCRV_remain_"$cluster_id
cp $outer_dir"asm_2nd_contigs/TCRV_remain_"$cluster_id"/contigs.fasta" $outer_dir"asm_2nd_contigs/TCRV_2nd_contigs_"$cluster_id".fasta"

# align alleles to 2nd_contigs
../bwa/bwa index     $outer_dir"asm_2nd_contigs/TCRV_2nd_contigs_"$cluster_id".fasta"
../bwa/bwa mem -t 16 $outer_dir"asm_2nd_contigs/TCRV_2nd_contigs_"$cluster_id".fasta" $outer_dir"asm_reads_and_alleles/TCRV_"$cluster_id"_allele.fasta" > $outer_dir"asm_2nd_alignment/align_TCRV_2nd_"$cluster_id".sam"

