import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def fasta_to_dict(fn_fasta):
    '''parse the fasta file into a dictionary'''
    # dict_name_SEQ {}
    #  - keys: seq_name
    #  - values: seq_SEQ
    dict_name_SEQ = {}
    with open(fn_fasta, 'r') as f_f:
        seq_name = ""
        seq_SEQ = ""
        for line in f_f:
            if line[0] == '>':
                if seq_name != "":
                    if dict_name_SEQ.get(seq_name):
                        print("WARNING! Duplicate sequence name:", seq_name)
                    dict_name_SEQ[seq_name] = seq_SEQ
                seq_name = line.strip()[1:].split(' ')[0]
                seq_SEQ = ""
            else:
                seq_SEQ += line.strip()
        if dict_name_SEQ.get(seq_name):
            print("WARNING! Duplicate sequence name:", seq_name)
        dict_name_SEQ[seq_name] = seq_SEQ
    return dict_name_SEQ


def parse_CIGAR(cigar):
    list_number = []
    list_operand = []
    tmp_num = 0
    for char in cigar:
        if char.isdigit():
            # since cigar will always starts in number
            tmp_num = tmp_num*10 + int(char)
        else:
            # since cigar will always ends in operator
            list_number.append(tmp_num)
            tmp_num = 0
            list_operand.append(char)
    return (list_number, list_operand)


def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def get_hamming_dist(
    str_a,
    str_b,
    thrsd=None,
    case_insensitive=True,
    try_both_orient=True
):
    '''
    Calculate the hamming distance between input strings

    Args:
    - str_a, str_b (Strings): equal-length strings
    - thrsd (Int/None):
        if None, return distance (Int)
        if Int, return True/False (distance is <= or > thrds)
    - case_insensitive (Bool): whether to consider cases
    - try_both_orient (Bool): set True to try both orientations

    Returns:
    Int or True/False
    '''

    assert (len(str_a) == len(str_b))
    if case_insensitive:
        str_a = str_a.upper()
        str_b = str_b.upper()

    dist = 0
    if thrsd == None:
        for i, a in enumerate(str_a):
            if a != str_b[i]:
                dist += 1

        if not try_both_orient:
            return dist
        else:
            r_dist = 0
            for i, a in enumerate(get_reverse_complement(str_a)):
                if a != str_b[i]:
                    r_dist += 1
            return min(dist, r_dist)
    else:
        dist = True
        for i, a in enumerate(str_a):
            if a != str_b[i]:
                dist += 1
            if dist > thrsd:
                dist = False
                break

        # if it's below thrsd, return True
        if dist != False:
            return True

        # if not, consider the other orientation
        if not try_both_orient:
            return False
        else:
            dist = True
            for i, a in enumerate(get_reverse_complement(str_a)):
                if a != str_b[i]:
                    dist += 1
                if dist > thrsd:
                    return False
            return True


def read_from_fastq(fn_fastq):
    '''
    Retrive information from a fastq file

    Inputs:
        - fn_fastq: input fastq file (.fq or .fastq)
    Outputs:
        - names: list of QNAMEs
        - seqs : list of SEQs
        - quals: list of QUALs
    '''
    ext = fn_fastq[fn_fastq.rfind('.') + 1:]
    assert ext in ['fastq', 'fq']

    with open(fn_fastq, 'r') as f:
        names = []
        seqs = []
        quals = []
        idx = 0
        for line in f:
            line = line.rstrip()
            if idx == 0:
                names.append(line[1:])
            elif idx == 1:
                seqs.append(line)
            elif idx == 2:
                quals.append(line)
            elif idx == 3:
                idx = 0
                continue
            idx += 1
    return names, seqs, quals

def read_from_file(fn_input):
    '''
    Reads sequences from a file.

    Inputs:
        - fn_input: input file name
            Supported formats:
            1. fasta, fa,
                typical fasta format
            2. txt
                each line represents a sequence
    Outputs:
        - list_seqs: list of sequences
    '''
    ext = fn_input[fn_input.rfind('.') + 1:]
    assert ext in ['txt', 'fasta', 'fa']

    list_seqs = []
    with open(fn_input, 'r') as f:
        for line in f:
            if ext in ['fasta', 'fa']:
                #: name
                if line[0] == '>':
                    continue
            line = line.rstrip().upper()
            
            #: doesn't allow nuc other than ACGT
            count_A = line.count('A')
            count_C = line.count('C')
            count_G = line.count('G')
            count_T = line.count('T')
            assert count_A + count_C + count_G + count_T == len(line)
            
            list_seqs.append(line)

    #: checks if all the seqs are equal in length
    for s in list_seqs:
        assert len(s) == len(list_seqs[0])
    
    return list_seqs

def read_allele_cluster(fn_cluster):
    '''
    Read the clustering result generated by clustal-omega
    
    Args:
    - fn_custer: file name of the clustering results

    Retuns:
    - dict_allele_cluster (Dict):
        keys are allele names, and values are cluster ids

    Example:
        From the following record
            1 Cluster 0: object 0 has index 17 (=seq M21626|TRAV14/DV4*01|Homo 290~len)▸   0000¬
        we have:
            {'TRAV14/DV4*01': '0'}
    '''
    f = open(fn_cluster, 'r')
    dict_allele_cluster = {}
    for line in f:
        cols = line.split()
        cluster_id = cols[1][:-1] # remove ":" in the last
        allele_name = cols[8].split('|')[1]
        dict_allele_cluster[allele_name] = cluster_id

    return dict_allele_cluster

def get_allele_cluster(dict_cluster, matched_allele_name):
    '''
    For an allele, get its cluster id
    '''
    if matched_allele_name.count('|'):
        return dict_cluster[matched_allele_name.split('|')[1]]
    return dict_cluster[matched_allele_name]

