import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fd', '--fn_data',
        help = 'allele data fasta'
    )
    parser.add_argument(
        '-fa', '--fn_annotated',
        help = 'annotated allele list'
    )
    parser.add_argument(
        '-o', '--fn_output',
        help = 'output fasta file with only annotated allele'
    )
    args = parser.parse_args()
    return args

def fetch_annotated(fn_data, fn_annotated):
    # gene = 'TRBV'
    
    list_annotated = []
    with open(fn_annotated, 'r') as f:
        for line in f:
            line = line.rstrip()
            allele = line.split(',')[0]
            list_annotated.append(allele)
    
    list_fetch = []
    write_flag = 0
    if fn_output:
        with open(fn_output, 'w') as f_o:
            with open(fn_data, 'r') as f:
                for line in f:
                    if line[0] == '>':
                        sline = line.strip()
                        allele = sline.split('|')[1]
                        list_fetch.append(allele)
                        if allele in list_annotated: # select the alleles with interest
                            f_o.write(line)
                            write_flag = 1
                        else:
                            write_flag = 0
                    elif write_flag: # output the sequence content with interest
                        f_o.write(line)
    else:
        print("\nWARNING: no output file specified!\n")

    set_annotated = set(list_annotated)
    set_fetched   = set(list_fetch)
    if len(set_annotated - set_fetched) > 0:
        print("\nWARNING: some names in the annotated list are not in the database!")
        for s in set_annotated - set_fetched:
            print(s)
        print()


if __name__ == '__main__':
    args = parse_args()
    fn_data = args.fn_data
    fn_annotated = args.fn_annotated
    fn_output = args.fn_output

    fetch_annotated(fn_data, fn_annotated)
