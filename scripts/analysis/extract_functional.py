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
    take in the read_depth report file
    and output only the functional genes
    """
    set_functional = set(list_functional)
    f = open(fn_report, 'r')
    for line in f:
        fields = line.strip().split()
        allele_name = fields[0]
        try:
            gene_name, allele_type = allele_name.split('*')
            if gene_name in set_functional:
                print(line.strip())
        except:
            print("Format Incorrect!", line)
    f.close()




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--read_depth_report', help='read depth report')
    parser.add_argument('-l', '--functional_list', help='functional gene list file')
    args = parser.parse_args()

    fn_report = args.read_depth_report
    fn_list   = args.functional_list

    list_functional = functional_list(fn_list)
    call_report(fn_report, list_functional)
