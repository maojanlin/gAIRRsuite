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
    parser.add_argument(
        '-o', '--fn_output',
        help = 'output file specifying the correctness of alleles [None]'
    )
    args = parser.parse_args()
    return args

'''
Read called and annotation files, filter by gene name, and store results in two sets.

Write analyzed results to output file if it is not empty.
'''
def read_files(fn_called, fn_annotated, gene, fn_output):
    list_blastn = []
    dict_blastn_allele_count = {}
    with open(fn_called, 'r') as f:
        for line in f:
            line = line.strip()
            allele = line.split('\t')[0]
            count = line.split('\t')[1]
            if line.count(gene) == 0:
                print ('Warning: gene {} not appearing in the blastn record: '.format(gene), line)
            else:
                list_blastn.append(allele)
                dict_blastn_allele_count[allele] = count

    list_annotated = []
    with open(fn_annotated, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.count(gene) == 0:
                print ('Warning: gene {} not appearing in the blastn annotation: '.format(gene), line)
            else:
                list_annotated.append(line)

    set_blastn = set(list_blastn)
    set_annotated = set(list_annotated)

    if fn_output:
        with open(fn_output, 'w') as f:
            
            for s in list_blastn:#.intersection(set_annotated):
                count = dict_blastn_allele_count[s]
                if s in set_annotated:
                    correctness = '3'
                else:
                    correctness = '2'
                f.write(s + '\t' + str(count) + '\t' + correctness + '\n')
            for s in set_annotated - set_blastn.intersection(set_annotated): #print those annotated but not found
                f.write(s + '\t' + '' + '\t' + '1' + '\n')

            '''
            for s in set_blastn.intersection(set_annotated):
                f.write(s + '\t' + '3' + '\n')
            set_in_blastn_not_annotated = set_blastn - set_blastn.intersection(set_annotated)
            for s in set_in_blastn_not_annotated:
                f.write(s + '\t' + '2' + '\n')
            set_not_in_blastn_but_annotated = set_annotated - set_blastn.intersection(set_annotated)
            for s in set_not_in_blastn_but_annotated:
                f.write(s + '\t' + '1' + '\n')
            '''
            
    return set_blastn, set_annotated, dict_blastn_allele_count

'''
Print results
'''
def print_results(set_blastn, set_annotated, dict_blastn_allele_count):
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
    fn_output = args.fn_output

    set_blastn, set_annotated, dict_blastn_allele_count = read_files(fn_called, fn_annotated, gene, fn_output)
    print_results(set_blastn, set_annotated, dict_blastn_allele_count)
