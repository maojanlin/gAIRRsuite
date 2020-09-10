import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fc', '--fn_contig',
        help = 'input contig fasta file'
    )
    parser.add_argument(
        '-id', '--cluster_id',
        help = 'the cluster id'
    )
    parser.add_argument(
        '-map', '--id_allele_map',
        help = 'cluster id to allele name mapping file'
    )
    parser.add_argument(
        '-fo', '--fo_fasta',
        help = 'output parsed contig fasta file'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_contig (fn_contig, fo_fasta, allele_name):
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
        print("NO CONTIG FILE", fn_contig, "FOUND!")
        

    if contig_name == "" or contig_SEQ == "":
        print("There is no contig assembled in", fn_contig, "!")
    else:
        f_o = open(fo_fasta, 'w')
        if allele_name:
            f_o.write('>' + allele_name + '_contig_' + contig_name.split('_')[5] + '\n')
        else:
            f_o.write(contig_name + '\n')
        f_o.write(contig_SEQ + '\n')
        f_o.close()
        print("Contig", fn_contig, "copy successfully!")



if __name__ == '__main__':
    args = parse_args()
    # input-output file
    fn_contig = args.fn_contig
    fo_fasta  = args.fo_fasta

    # parameter for rename the contigs
    cluster_id = args.cluster_id
    id_allele_map = args.id_allele_map

    allele_name = None
    if cluster_id and id_allele_map:
        f_map = open(id_allele_map, 'r')
        for line in f_map:
            allele_id, allele_name = (line.strip()).split(',')
            if allele_id == cluster_id:
                break
        f_map.close()
    parse_contig(fn_contig, fo_fasta, allele_name)
