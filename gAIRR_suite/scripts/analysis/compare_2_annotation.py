import argparse


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




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a1', '--annotation_report_1', help='the 1st annotation file comparing')
    parser.add_argument('-a2', '--annotation_report_2', help='the 2nd annotation file comparing')
    args = parser.parse_args()

    fn_r1 = args.annotation_report_1
    fn_r2 = args.annotation_report_2

    list_known_1, list_novel_1 = parse_annotation_report(fn_r1)
    list_known_2, list_novel_2 = parse_annotation_report(fn_r2)
    set_known_1 = set(list_known_1)
    set_known_2 = set(list_known_2)

    print("For known alleles:")
    print("There are", len(set_known_1.intersection(set_known_2)), "alleles in common")
    print(fn_r1, "only alleles:\n", sorted(set_known_1 - set_known_2))
    print(fn_r2, "only alleles:\n", sorted(set_known_2 - set_known_1))
    print("Novel alleles in", fn_r1, ":\n", sorted(list_novel_1))
    print("Novel alleles in", fn_r2, ":\n", sorted(list_novel_2))
