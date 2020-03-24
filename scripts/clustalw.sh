#use clustalw to msa annotated_allele.fasta
./clustalw2 ./experiment/TRBV_filtered.fasta
# should produce alignment .aln file and tree .dnd file in the target directory

# parse .dnd file
python3 ../immunogenomics/script/parse_tree.py -ft ./experiment/TRBV_filtered.dnd -o ./experiment/TRBV_parsed_tree.dnd

# plot tree
Rscript ../immunogenomics/script/plot_tree.r -t ./experiment/TRBV_parsed_tree.dnd -a ./experiment/trbv_top80.labelled.txt -s 1 -o ./experiment/tree.pdf
