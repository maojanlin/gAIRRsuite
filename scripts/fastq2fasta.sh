# arg 1: input filename (fastq we'd like to convert)
# arg 2: output filename (converted fasta)
sed -n '1~4s/^@/>/p;2~4p' $1 > $2
