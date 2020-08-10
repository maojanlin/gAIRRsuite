import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fc', '--fn_contig',
        help = 'input contig fasta file'
    )
    parser.add_argument(
        '-fo', '--fo_fasta',
        help = 'output parsed contig fasta file'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_contig (fn_contig, fo_fasta):
    contig_name = ""
    contig_SEQ  = ""
    try:
        with open(fn_contig, 'r') as f_o:
            for line in f_o:
                if line[0] == '>':
                    if contig_name == "":
                        contig_name = line.strip()
                    else:
                        break
                else:
                    contig_SEQ += line.strip()
    except FileNotFoundError:
        print("NO CONTIG FILE", fn_contig, "FOUND!!!")
        

    if contig_name == "" or contig_SEQ == "":
        print("There is no contig in", fn_contig, "!!!")
    else:
        f_o = open(fo_fasta, 'w')
        f_o.write(contig_name + '\n')
        f_o.write(contig_SEQ + '\n')
        f_o.close()
        print("Contig", fn_contig, "copy successfully!")



if __name__ == '__main__':
    args = parse_args()
    fn_contig = args.fn_contig
    fo_fasta  = args.fo_fasta

    parse_contig(fn_contig, fo_fasta)
