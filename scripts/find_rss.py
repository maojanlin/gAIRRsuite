fn = 'read-RA_si-AGCTATCA_lane-001-chunk-002.fastq'
with open(fn, 'r') as f:
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
cnt_hep = 0
cnt_non = 0
cnt_both = 0
f_out = open('exact_pm2_forward.names', 'w')
for i, s in enumerate(seqs):
    has_hep = s.count('CACAGTG')
    has_non = s.count('ACAAAAACC')
    if has_hep > 0:
        cnt_hep += 1
    if has_non > 0:
        cnt_non += 1
    if has_hep and has_non:
        cnt_both += 1
        pos_hep = s.find('CACAGTG')
        pos_non = s.find('ACAAAAACC')
        # print (pos_hep, pos_non, pos_non - pos_hep)
        #: 12/23
        if pos_non - pos_hep in [17,18,19,20,21, 28,29,30,31,32]:
            f_out.write(names[i] + '\n')
            print ('inrange', names[i])