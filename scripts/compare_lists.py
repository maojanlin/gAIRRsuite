import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fc', '--fn_called',
        help = 'called allele list'
    )
    parser.add_argument(
        '-fa', '--fn_annotated',
        help = 'annotated allele list'
    )
    parser.add_argument(
        '-g', '--gene',
        help = 'gene name of interest'
    )
    args = parser.parse_args()
    return args

def compare_files(fn_called, fn_annotated, gene):
    # gene = 'TRBV'
    list_blastn = []
    dict_blastn_allele_count = {}
    with open(fn_called, 'r') as f:
        for line in f:
            line = line.strip()
            allele = line.split(':')[0]
            count = line.split(':')[1]
            if line.count(gene) == 0:
                print ('warning:', line)
            else:
                list_blastn.append(allele)
                dict_blastn_allele_count[allele] = count

    list_annotated = []
    with open(fn_annotated, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.count(gene) == 0:
                print ('warning:', line)
            else:
                list_annotated.append(line)

    set_blastn = set(list_blastn)
    set_annotated = set(list_annotated)

    print ()
    print ('Size of blastn alleles:', len(set_blastn))
    print ('Size of annotated alleles:', len(set_annotated))
    print ('Size of intersection:', len(set_blastn.intersection(set_annotated)))
    print ('Size of union:', len(set_blastn.union(set_annotated)))

    set_in_blastn_not_annotated = set_blastn - set_blastn.intersection(set_annotated)
    print ()
    print ('Found by blastn, not annotated:')
    for s in set_in_blastn_not_annotated:
        print (s, dict_blastn_allele_count[s])
    set_not_in_blastn_but_annotated = set_annotated - set_blastn.intersection(set_annotated)
    print ()
    print ('Not found by blastn, but annotated:')
    for s in set_not_in_blastn_but_annotated:
        print (s)


if __name__ == '__main__':
    args = parse_args()
    fn_called = args.fn_called
    fn_annotated = args.fn_annotated
    gene = args.gene

    compare_files(fn_called, fn_annotated, gene)
