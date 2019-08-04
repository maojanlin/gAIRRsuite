#: arg1: input file (FASTQ)
#: arg2: output file (FASTQ)
bgzip -cd $1 | grep -A 2 -B 1 --no-group-separator "CACAGTG\|ACAAAAACC\|GGTTTTTGT" > $2
