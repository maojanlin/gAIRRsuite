# make a valid name from real allele name
echo "Allele_name: "$1
valid_name=$(echo $1 | sed 's/*/_/g')
valid_name=$(echo $valid_name | sed 's/\//_/g')
valid_name=$(echo $valid_name | sed 's/(/_/g')
valid_name=$(echo $valid_name | sed 's/)/_/g')

# parse the sam file to find out the H1 contig and clip region of the allele
sam_output=$(python3 ../../immunogenomics/scripts/fetch_sam.py -fs ./support_reads_analysis/IG_alleles_H1.sam -fa $1)
contig_name=$(echo $sam_output | cut -d" " -f1)
clip_region=$(echo $sam_output | cut -d" " -f2)
allele_SEQ=$(echo $sam_output | cut -d" " -f3)
echo "In NA12878-H1:"
echo "Contig_name: "$contig_name", clip_region: "$clip_region

# fetch the clip region reads from align_to_assembly bam file
samtools view ../../../align_to_asm/20200514_NA12878_aligntoasm_H1/NA12878_align_to_H1asm.sorted.bam $contig_name":"$clip_region > "./support_reads_analysis/"$valid_name"_interval_H1.txt"

# parse the genome file to get the contig sequence
python3 ../../immunogenomics/scripts/fetch_chromosome.py -fg ../../../asm/NA12878/NA12878-H1.fa -fc $contig_name -fo "./support_reads_analysis/"$valid_name"_contig_H1.fasta"

# call where allele and contig disagree or site with SNP
python3 ../../immunogenomics/scripts/parse_sam_haplotyping.py -fs "./support_reads_analysis/"$valid_name"_interval_H1.txt" --interval $clip_region -fas $allele_SEQ -fc "./support_reads_analysis/"$valid_name"_contig_H1.fasta"


# parse the sam file to find out the H2 contig and clip region of the allele
sam_output=$(python3 ../../immunogenomics/scripts/fetch_sam.py -fs ./support_reads_analysis/IG_alleles_H2.sam -fa $1)
contig_name=$(echo $sam_output | cut -d" " -f1)
clip_region=$(echo $sam_output | cut -d" " -f2)
allele_SEQ=$(echo $sam_output | cut -d" " -f3)
echo "In NA12878-H2:"
echo "Contig_name: "$contig_name", clip_region: "$clip_region
# fetch the clip region reads from align_to_assembly bam file
samtools view ../../../align_to_asm/20200514_NA12878_aligntoasm_H2/NA12878_align_to_H2asm.sorted.bam $contig_name":"$clip_region > "./support_reads_analysis/"$valid_name"_interval_H2.txt"
# parse the genome file to get the contig sequence
python3 ../../immunogenomics/scripts/fetch_chromosome.py -fg ../../../asm/NA12878/NA12878-H2.fa -fc $contig_name -fo "./support_reads_analysis/"$valid_name"_contig_H2.fasta"
# call where allele and contig disagree or site with SNP
python3 ../../immunogenomics/scripts/parse_sam_haplotyping.py -fs "./support_reads_analysis/"$valid_name"_interval_H2.txt" --interval $clip_region -fas $allele_SEQ -fc "./support_reads_analysis/"$valid_name"_contig_H2.fasta"
