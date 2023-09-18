import argparse


def functional_list(fn_list):
    """
    extract the list log and return a list containing all the gene names
    """
    list_functional = []
    f = open(fn_list, 'r')
    for line in f:
        list_functional.append(line.strip())
    f.close()
    return list_functional


def call_report(
        fn_report       :str,
        list_functional :list
    ) -> None:
    """
    take in the annotation report file
    and output only the functional genes
    """
    dict_functional = {}
    f = open(fn_report, 'r')
    for line in f:
        fields = line.strip().split(',')
        if len(fields) == 1:
            continue
        allele_name = fields[0]
        try:
            gene_name, allele_type = allele_name.split('*')
            if len(fields) > 4:
                allele_name += '_novel'
            if dict_functional.get(gene_name):
                dict_functional[gene_name].add(allele_name)
            else:
                dict_functional[gene_name] = {allele_name}
        except:
            print("Format Incorrect!", line)
    f.close()

    if list_functional:
        for gene_name in list_functional:
            if dict_functional.get(gene_name):
                list_allele = dict_functional[gene_name]
                for name in sorted(list_allele):
                    print(name)
    else:
        for gene_name, list_allele in sorted(dict_functional.items()):
            for name in sorted(list_allele):
                print(name)





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--annotation_report', help='gAIRR-annotate report')
    parser.add_argument('-l', '--functional_list', help='functional gene list file')
    args = parser.parse_args()

    fn_report = args.annotation_report
    fn_list   = args.functional_list

    if fn_list:
        list_functional = functional_list(fn_list)
    else:
        list_functional = None
    call_report(fn_report, list_functional)

