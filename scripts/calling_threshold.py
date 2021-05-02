import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-dp', '--fn_depth_report',
        help = 'input read-depth calling report'
    )
    args = parser.parse_args()
    return args


def find_thresh(list_depth):
    p_depth = 1
    thresh = 1
    for idx, (allele_name, depth) in enumerate(list_depth):
        if 'TRAV8-5*01' in allele_name: # the special allele
            continue
        if depth / p_depth < 0.75:
            thresh = depth + 1
            print("THreshold!!", thresh)
            return thresh
        p_depth = depth
    return thresh


def thresh_divide(list_depth, thresh):
    total_num = 0
    novel_num = 0
    flag = True
    p_depth = 1
    for allele_name, depth in list_depth:
        if depth < thresh:
            if flag:
                flag = False
                print("------------- thresh ----------------")
            pass
        else:
            total_num += 1
            if 'novel' in allele_name:
                novel_num += 1
        print(allele_name, depth)#, depth/p_depth, sep='\t\t')
        p_depth = depth
        if depth == 0:
            return total_num, novel_num
    return total_num, novel_num



if __name__ == '__main__':
    args = parse_args()
    fn_depth_report = args.fn_depth_report

    f_n = open(fn_depth_report, 'r')
    list_depth = [] # list_depth = [(name1,depth1), (name2,depth2), ... ]
    for line in f_n:
        fields = line.split()
        allele_name = fields[0]
        depth = float(fields[1])
        list_depth.append((allele_name, depth))
    f_n.close()

    thresh = find_thresh(list_depth)
    total_num, novel_num = thresh_divide(list_depth, thresh)

    print("\n========= Summary ===========")
    print("Total AIRRCall alleles:", total_num)
    print("Novel AIRRCall alleles:", novel_num)
    
