#use clustalw to msa annotated_allele.fasta
./clustalw2 ./experiment/TRBV_filtered.fasta
# should produce alignment .aln file and tree .dnd file in the target directory

# parse .dnd file
python3 ../script/parse_tree.py -ft ./experiment/TRBV_filtered.dnd -o ./experiment/TRBV_parsed_tree.dnd
