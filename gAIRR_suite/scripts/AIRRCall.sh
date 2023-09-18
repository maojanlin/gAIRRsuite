workspace="target_call"
path_SPAdes="spades.py"
list_allele_name="TCRJ TCRV BCRV BCRJ"
allele_dir="./example/material/"
allele_suffix="_alleles_parsed.fasta"
person_name="HG002-part"
read_path_1="./example/samples/HG002_part_gAIRR-seq_R1.fastq.gz"
read_path_2="./example/samples/HG002_part_gAIRR-seq_R2.fastq.gz"
thread=16

# environment settings
mkdir -p ${workspace}
mkdir -p ${workspace}/${person_name}
workspace=${workspace}/${person_name}

# call novel alleles
for allele_name in ${list_allele_name}; do
    allele_path=${allele_dir}${allele_name}${allele_suffix}
    echo "[AIRRCall] ${person_name} ${allele_name} novel allele calling..."
    ./scripts/novel_allele.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2} ${thread}
    
    allele_path=${workspace}/${person_name}_${allele_name}_novel/${allele_name}_with_novel.fasta
    echo "[AIRRCall] ${person_name} ${allele_name} allele calling..."
    ./scripts/allele_calling.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2} ${thread}
    
    echo "[AIRRCall] ${person_name} ${allele_name} flanking sequence calling..."
    ./scripts/flanking_sequence.sh ${workspace} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2} ${path_SPAdes} ${thread}
done
