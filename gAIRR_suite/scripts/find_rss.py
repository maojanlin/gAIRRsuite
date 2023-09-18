import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--fn_fastq',
        help='filtered fastq file'
    )
    parser.add_argument(
        '-o', '--fn_out',
        help='output file (names)'
    )
    parser.add_argument(
        '-r', '--tol_range', type=int, default=0,
        help='tolerance of distance between hep and non [1]'
    )
    args = parser.parse_args()
    return args

#fn = 'read-RA_si-AGCTATCA_lane-001-chunk-002-rss_filtered.fastq'
#f_out = open(fn + '-rss.names', 'w')

def write_when_find_both(
    seq,
    idx,
    first,
    second,
    names,
    f_out,
    tol_range
):
    cnt_both = 0
    has_first = seq.count(first)
    has_second = seq.count(second)
    if has_first and has_second:
        cnt_both += 1
        pos_first = seq.find(first)
        pos_second = seq.find(second)
        # print (pos_first, pos_second, pos_second - pos_first)
        #: 12/23
        diff_and_range = \
                list(
                    range(12 + len(first) - tol_range, 12 + len(first) + tol_range)
                ) + \
                list(
                    range(23 + len(first) - tol_range, 23 + len(first) + tol_range)
                )
        #if pos_second - pos_first in [17,18,19,20,21, 28,29,30,31,32]:
        if pos_second - pos_first in diff_and_range:
            f_out.write(names[idx] + '\n')
            print ('inrange', names[idx])
    return cnt_both

if __name__ == '__main__':
    args = parse_args()
    fn_fastq = args.fn_fastq
    fn_out = args.fn_out
    f_out = open(fn_out, 'w')
    tol_range = args.tol_range

    with open(fn_fastq, 'r') as f:
        names = []
        seqs = []
        idx = 0
        for line in f:
            line = line.rstrip()
            if idx == 0:
                names.append(line[1:])
            elif idx == 1:
                seqs.append(line)
            elif idx == 3:
                idx = 0
                continue
            idx += 1
    
    cnt_forward = 0
    cnt_reverse = 0
    
    for i, s in enumerate(seqs):
        seq_hep = 'CACAGTG'
        seq_non = 'ACAAAAACC'
        #: forward
        cnt_forward += write_when_find_both(
            seq=s,
            idx=i,
            first=seq_hep,
            second=seq_non,
            names=names,
            f_out=f_out,
            tol_range=tol_range
        )
        #: reverse
        seq_hep = 'CACAGTG'
        seq_non = 'GGTTTTTGT'
        cnt_reverse += write_when_find_both(
            seq=s,
            idx=i,
            first=seq_non,
            second=seq_hep,
            names=names,
            f_out=f_out,
            tol_range=tol_range
        )
    
    print ('forward', cnt_forward)
    print ('reverse', cnt_reverse)
