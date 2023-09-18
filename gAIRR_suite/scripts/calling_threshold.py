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
    """
    Threshold algorithm:
        - threshold should <= average # considering 0s are filtered
        - a sliding window size 3
    """
    list_d = [pair_info[1] for pair_info in list_depth if 'TRAV8-5*01' not in pair_info[0]] # rebuild a list without TRAV8-5*01 (outlier) stuff 
    avg_d  = sum(list_d)/len(list_d)
    list_d = [list_d[0]]*3 + list_d + [list_d[-1]]*2                                        # padding max value and zeros
    #print(list_d)
    list_value = []     # the division value between two windows
    for idx in range(3,len(list_d)-2):
        window   = (list_d[idx] + list_d[idx+1]*0.25 + list_d[idx+2]*0.1)
        p_window = (list_d[idx-3]*0.1 + list_d[idx-2]*0.25 + list_d[idx-1])
        list_value.append((((window+2)/(p_window+2))*((window+2)/(window+0.5)),idx))
        #print(idx-3, list_d[idx], format(((window+2)/(p_window+2))*((window+2)/(window+0.5)), '.3f'))
    #print(sorted(list_value))
    sorted_value = sorted(list_value)
    if len(sorted_value) <= 1:
        return 1
    thresh_id = -1
    if sorted_value[0][0]*2 < sorted_value[1][0]: # absolute winner
        thresh_id = sorted_value[0][1]
    else:                                         # if there are similar candidate, use the old method
        for idx in range(3,len(list_d)-2):        # old method: find single site has slope greater than 0.73
            if list_d[idx] / list_d[idx-1] < 0.73 and list_d[idx] <= avg_d:
                thresh_id = idx
                break
    thresh = list_d[thresh_id] + 1
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
                print("------------- thresh: " + str(thresh) + " ----------------")
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
        if depth == 0: # keep the first 0
            break
    f_n.close()

    thresh = find_thresh(list_depth)
    total_num, novel_num = thresh_divide(list_depth, thresh)

    print("\n========= Summary ===========")
    print("Total AIRRCall alleles:", total_num)
    print("Novel AIRRCall alleles:", novel_num)
    
