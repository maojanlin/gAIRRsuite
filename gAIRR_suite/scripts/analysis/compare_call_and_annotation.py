import argparse
import sys



def parse_annotation_report(fn_report):
    list_known = []
    list_novel = []
    f = open(fn_report, 'r')
    for line in f:
        allele_name = line.strip()
        if 'novel' in allele_name:
            list_novel.append(allele_name)
        else:
            list_known.append(allele_name)
    f.close()
    return list_known, list_novel


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def check_call(fn_call, list_annotate_known, list_annotate_novel):
    set_annotate_known = set(list_annotate_known)
    list_call_novel = []
    f = open(fn_call, 'r')
    for line in f:
        fields = line.strip().split()
        allele_name = fields[0]
        if 'novel' in allele_name:
            if float(fields[1]) > 0:
                print(line.strip() + '\t' + '?')
                list_call_novel.append(allele_name)
            else:
                print(line.strip() + '\t' + '2')
        elif allele_name in set_annotate_known:
            print(line.strip() + '\t' + '3')
        else:
            print(line.strip() + '\t' + '2')
    f.close()

    eprint("Novel alleles in annotation:")
    eprint(list_annotate_novel)
    eprint("Novel alleles in calling:")
    eprint(list_call_novel)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-an', '--annotation_report', help='the annotation file comparing')
    parser.add_argument('-rc', '--read_depth_call', help='the read depth file comparing')
    args = parser.parse_args()

    fn_annotate = args.annotation_report
    fn_call     = args.read_depth_call

    list_known, list_novel = parse_annotation_report(fn_annotate)
    check_call(fn_call, list_known, list_novel)
    

