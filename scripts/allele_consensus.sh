#!/bin/bash
# the out most directory
outer_dir="HGSVC_annotation/"
list_allele_name="TCRV TCRD_plusHep TCRJ BCRV BCRD_plusHep BCRJ"
name_log="consensus_name_HGSVC.log"

echo "[genAIRR] Building the Consensus Database in ${outer_dir}"
for allele_name in ${list_allele_name}; do
    python3 scripts/allele_consensus.py -fn  ${name_log} \
                                        -fr  ${outer_dir}database_novel_${allele_name}.tsv \
                                        -ff  ${outer_dir}database_novel_${allele_name}.fasta \
                                        -for ${outer_dir}database_novel_${allele_name}_consensus.tsv \
                                        -fof ${outer_dir}database_novel_${allele_name}_consensus.fasta
    
    python3 scripts/allele_consensus.py -fn  ${name_log} \
                                        -fr  ${outer_dir}database_flanking_${allele_name}.tsv \
                                        -ff  ${outer_dir}database_flanking_${allele_name}.fasta \
                                        -for ${outer_dir}database_flanking_${allele_name}_consensus.tsv \
                                        -fof ${outer_dir}database_flanking_${allele_name}_consensus.fasta
done
echo "[genAIRR] Consensus Database Finished!"
