workspace="target_call"
list_allele_name="TCRJ TCRV BCRV"
allele_dir="../genomeData/"
allele_suffix="_alleles_parsed.fasta"
person_name="NA12878"
read_path_1="../captured_reads/NA12878_S46_L001_R1_001.fasta"
read_path_2="../captured_reads/NA12878_S46_L001_R2_001.fasta"

# environment settings
mkdir -p ${workspace}

# call novel alleles
for allele_name in ${list_allele_name}; do
    allele_path=${allele_dir}${allele_name}${allele_suffix}
    echo "[AIRRCall] ${person_name} ${allele_name} novel allele calling..."
    ./scripts/novel_allele.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
    
    allele_path=${workspace}/${person_name}_${allele_name}_novel/${allele_name}_with_novel.fasta
    echo "[AIRRCall] ${person_name} ${allele_name} allele calling..."
    ./scripts/allele_calling.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
    
    echo "[AIRRCall] ${person_name} ${allele_name} flanking sequence calling..."
    ./scripts/flanking_sequence.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
done
