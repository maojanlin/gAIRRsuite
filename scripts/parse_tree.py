import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-ft', '--fn_tree',
        help = 'dnd (or Newick format) tree file'
    )
    parser.add_argument(
        '-o', '--fn_output',
        help = 'output tree file with the most concise name'
    )
    args = parser.parse_args()
    return args

def parse_tree(fn_tree):
    
    list_annotated = []
    with open(fn_output, 'w') as f_o:
        with open(fn_tree, 'r') as f:
            for line in f:
                line = line.rstrip()
                if len(line.split('|')) > 2:
                    allele = line.split(',')[0]
                    parsed_string = line.split('|')[1] + ':' + line.split(':')[1] + '\n'
                    f_o.write(parsed_string)
                else:
                    f_o.write(line + '\n')



if __name__ == '__main__':
    args = parse_args()
    fn_tree = args.fn_tree
    fn_output = args.fn_output

    parse_tree(fn_tree)
