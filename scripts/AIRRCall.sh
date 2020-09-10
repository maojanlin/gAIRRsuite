#allele_path="../plot_tree/TCRJ_alleles.fasta"
#read_path_1="CHM13_filtered_R1.fasta"
#read_path_2="CHM13_filtered_R2.fasta"

workspace="target_call"
allele_path="../plot_tree/TCRJ_alleles.fasta"
allele_name="TCRJ"
person_name="NA12878"
read_path_1="NA12878_S46_R1.fasta"
read_path_2="NA12878_S46_R2.fasta"

# environment settings
mkdir -p ${workspace}

# call novel alleles
echo "[AIRRCall] ${allele_name} novel allele calling..."
#./novel_allele.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
allele_path=${workspace}/${person_name}_${allele_name}_novel/${allele_name}_with_novel.fasta

echo "[AIRRCall] ${allele_name} allele calling..."
#./allele_calling.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}

echo "[AIRRCall] ${allele_name} flanking sequence calling..."
./flanking_sequence.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
