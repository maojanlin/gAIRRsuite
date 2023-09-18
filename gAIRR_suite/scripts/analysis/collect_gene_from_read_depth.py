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
        fn_read_depth    :str,
        set_functional :list
    ) -> None:
    """
    take in the annotation report file
    and output only the functional genes
    """
    dict_functional = {}
    for gene_name in set_functional:
        dict_functional[gene_name] = []

    # take in read depth report allele types
    f = open(fn_read_depth, 'r')
    for line in f:
        fields = line.strip().split()
        if 'thresh' in line:
            break
        allele_name = fields[0]
        gene_name, allele_type = allele_name.split('*')
        if dict_functional.get(gene_name) != None:
            dict_functional[gene_name].append(allele_type)
    f.close()
    return dict_functional



def print_report(dict_genes, dict_functional):
    for gene, list_functional in dict_genes.items():
        print("# " + gene)
        for gene_name in list_functional:
            allele_type = dict_functional[gene_name]
            print(gene_name, allele_type)






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-rd', '--read_depth_report', help='gAIRR-call read depth report')
    parser.add_argument('-fl', '--functional_list', nargs='+', default=[], help='functional gene list file')
    args = parser.parse_args()

    fn_read_depth = args.read_depth_report
    fn_list       = args.functional_list

    dict_genes, set_functional = functional_list(fn_list)
    dict_functional = call_report(fn_read_depth, set_functional)
    print_report(dict_genes, dict_functional)

