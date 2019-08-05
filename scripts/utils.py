def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

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