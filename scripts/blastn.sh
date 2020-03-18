# Build BLASTn database
ncbi-blast-2.10.0+/bin/makeblastdb -in ../../20200317_filtered_V_alleles_for_probe_design/TRBV_filtered.fasta -dbtype nucl
# Convert FASTQ to FASTA
sh ../../fastq2fasta.sh NA12878_S46_L001_R1_001.fastq NA12878_S46_L001_R1_001.fasta
# Run BLASTn
ncbi-blast-2.10.0+/bin/blastn -query experiments/NA12878_S46_L001_R1_001.fasta -db ../../20200317_filtered_V_alleles_for_probe_design/TRBV_filtered.fasta -task blastn -num_threads 8 -outfmt "7 qseqid sseqid evalue bitscore" > experiments/NA12878_S46_L001_R1_001.blastn.out.txt
