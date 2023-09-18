import argparse


def functional_list(fn_list):
    """
    extract the list log and return a list containing all the gene names
    """
    dict_genes = {}
    set_functional  = set()
    for fn_name in fn_list:
        target_gene = fn_name.split('.')[0]
        list_functional = []
        f = open(fn_name, 'r')
        for line in f:
            gene_name = line.strip()
            list_functional.append(gene_name)
            set_functional.add(gene_name)
        f.close()
        dict_genes[target_gene] = list_functional
    return dict_genes, set_functional


def call_report(
        fn_bed_1        :str,
        fn_bed_2        :str,
        set_functional  :list
    ) -> None:
    """
    take in the annotation report file
    and output only the functional genes
    """
    dict_functional = {}
    for gene_name in set_functional:
        dict_functional[gene_name] = [[],[]]

    # take in haplotype 1 allele types
    f = open(fn_bed_1, 'r')
    for line in f:
        fields = line.strip().split()
        if len(fields) < 4:
            continue
        allele_name = fields[3]
        list_allele = allele_name.split(';')
        for allele in list_allele:
            gene_name, allele_type = allele.split('*')
            if dict_functional.get(gene_name):
                dict_functional[gene_name][0].append(allele_type)
    f.close()

    # take in haplotype 2 allele types
    f = open(fn_bed_2, 'r')
    for line in f:
        fields = line.strip().split()
        if len(fields) < 4:
            continue
        allele_name = fields[3]
        list_allele = allele_name.split(';')
        for allele in list_allele:
            gene_name, allele_type = allele.split('*')
            if dict_functional.get(gene_name):
                dict_functional[gene_name][1].append(allele_type)
    f.close()

    return dict_functional


def print_report(dict_genes, dict_functional):
    for gene, list_functional in dict_genes.items():
        print("# " + gene)
        for gene_name in list_functional:
            allele_H1, allele_H2 = dict_functional[gene_name]
            print(gene_name, allele_H1, allele_H2)



def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-b1', '--annotation_bed_1', help='gAIRR-annotate H1 bed file')
    parser.add_argument('-b2', '--annotation_bed_2', help='gAIRR-annotate H2 bed file')
    parser.add_argument('-fl', '--functional_list', nargs='+', default=[], help='functional gene list file')
    args = parser.parse_args(arguments)

    fn_bed_1 = args.annotation_bed_1
    fn_bed_2 = args.annotation_bed_2
    fn_list  = args.functional_list

    dict_genes, set_functional = functional_list(fn_list)
    dict_functional = call_report(fn_bed_1, fn_bed_2, set_functional)
    print_report(dict_genes, dict_functional)





if __name__ == "__main__":
    main()
