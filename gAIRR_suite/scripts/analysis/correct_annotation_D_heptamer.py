import pysam
import argparse


def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-in', '--input_bed', help='input bed file for D gene correction')
    #parser.add_argument('-out', '--output', help='the output corrected bed file')
    args = parser.parse_args(arguments)
    
    fn_in   = args.input_bed
    #fn_out  = args.output

    with open(fn_in) as f:
        for line in f:
            fields = line.split()
            if len(fields) < 4:
                print(line.strip())
                continue
            gene_locus = fields[3][:4]
            if gene_locus == 'IGHD' or gene_locus == 'TRDD' or gene_locus == 'TRBD':
                start_pos = int(fields[1])
                end_pos   = int(fields[2])
                delimiter = line[len(fields[0])]
                print(delimiter.join([fields[0], str(start_pos+7), str(end_pos-7), fields[3].strip()]))
            else:
                print(line.strip())



if __name__ == "__main__":
    main()

