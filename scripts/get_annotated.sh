### Script to process annotations, see below example:
# cut -d ',' -f 1 ../../../20200317_NA12878_TCRV_annotated_alleles/TRBV_align_NA12878_asm_H* > trbv_all.txt
# Usage get_annotated.sh <alleles.csv> <output.txt>
cut -d ',' -f 1 $1 > $2
