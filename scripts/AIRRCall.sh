#allele_path="../plot_tree/TCRJ_alleles.fasta"
#read_path_1="CHM13_filtered_R1.fasta"
#read_path_2="CHM13_filtered_R2.fasta"

workspace="target_call"
allele_path="../genomeData/BCRV_alleles_parsed.fasta"
allele_name="BCRV"
person_name="NA12878"
read_path_1="../../naechyun/blast/experiments/NA12878_S46_L001_R1_001.fasta"
read_path_2="../../naechyun/blast/experiments/NA12878_S46_L001_R2_001.fasta"

# environment settings
mkdir -p ${workspace}

# call novel alleles
echo "[AIRRCall] ${person_name} ${allele_name} novel allele calling..."
./scripts/novel_allele.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
allele_path=${workspace}/${person_name}_${allele_name}_novel/${allele_name}_with_novel.fasta

echo "[AIRRCall] ${person_name} ${allele_name} allele calling..."
./scripts/allele_calling.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}

echo "[AIRRCall] ${person_name} ${allele_name} flanking sequence calling..."
./scripts/flanking_sequence.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
