import argparse
from utils import fasta_to_dict

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fn', '--fn_name_list',
        help = 'input name stating who has multiple assemblies'
    )
    parser.add_argument(
        '-fr', '--fn_report',
        help = 'input allele report file'
    )
    parser.add_argument(
        '-ff', '--fn_fasta',
        help = 'input allele fasta file'
    )
    parser.add_argument(
        '-for', '--fo_consensus_report',
        help = 'output report with consensus allele calling'
    )
    parser.add_argument(
        '-fof', '--fo_consensus_fasta',
        help = 'output fasta file of consensus allele sequence'
    )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    fn_name_list = args.fn_name_list
    fn_report = args.fn_report
    fn_fasta  = args.fn_fasta
    fo_consensus_report = args.fo_consensus_report
    fo_consensus_fasta  = args.fo_consensus_fasta

    # dict_sample_to_consensus{}
    # - keys: sample_name
    # - value: consensus_name
    dict_sample_to_consensus = {}
    # dict_sample_num{}
    # - keys: consensus_name
    # - values: number of samples
    dict_sample_num = {}
    
    f_n = open(fn_name_list, 'r')
    for line in f_n:
        fields = line.split()
        consensus_name = fields[0]
        list_sample_name = fields[1:]
        
        dict_sample_num[consensus_name] = len(list_sample_name)
        for sample_name in list_sample_name:
            dict_sample_to_consensus[sample_name] = consensus_name
    f_n.close()

    #print(dict_sample_to_consensus)
    #print(dict_sample_num)

    # dict_report{}
    # - keys: allele_name
    # - values: [sample_1, sample_2, sample_3, ...]
    dict_report = {}
    dict_fasta = fasta_to_dict(fn_fasta)
    f_r = open(fn_report, 'r')
    f_r.readline()
    for line in f_r:
        fields = line.split()
        allele_name = fields[0]
        list_sample = fields[2].split(',')
        for idx, ele in enumerate(list_sample):
            new_sample_name = ele.split('_')[0]
            list_sample[idx] = new_sample_name
        dict_report[allele_name] = list_sample

    # files are prepared
    for allele, list_sample in sorted(dict_report.items()):
        dict_num = {}
        for idx, sample in reversed(list(enumerate(list_sample))):
            if dict_sample_to_consensus.get(sample):
                list_sample.pop(idx)
                consensus_name = dict_sample_to_consensus[sample]
                if dict_num.get(consensus_name):
                    dict_num[consensus_name] += 1
                else:
                    dict_num[consensus_name] = 1

        for consensus_name, num in dict_num.items():
            if num > dict_sample_num[consensus_name]/2: # more than half consent
                list_sample.append(consensus_name)
        # check if there are still sample support the allele
        if len(list_sample) == 0:
            dict_report.pop(allele)
            dict_fasta.pop(allele)

    # polish allele name
    dict_old_name_2_new = {}
    old_name = None
    flank_idx = -1
    novel_idx = -1
    old_n_id = -1
    for allele_name in sorted(dict_report.keys()):
        if '/n' in allele_name and '/f' in allele_name: # complicated /n01/f01 cases
            current_name = allele_name[:-6]
            current_n_id = allele_name[-6:-4]
            if current_name == old_name: # same gene
                if old_n_id == current_n_id: # same novel allele
                    flank_idx += 1
                else: # different novel allele
                    flank_idx = 1
                    novel_idx += 1
            else: # change gene
                flank_idx = 1
                novel_idx = 1
                old_name = current_name
                old_n_id = current_n_id
            new_name = current_name + str(novel_idx).zfill(2) + '/f' + str(flank_idx).zfill(2)
            dict_old_name_2_new[allele_name] = new_name
        else: # /n01 only or /f01 only cases
            current_name = allele_name[:-2]
            if current_name == old_name:
                flank_idx += 1
            else:
                flank_idx = 1
                old_name = current_name
            new_name = current_name + str(flank_idx).zfill(2)
            dict_old_name_2_new[allele_name] = new_name


    f_of = open(fo_consensus_fasta,  'w')
    f_or = open(fo_consensus_report, 'w')
    f_or.write('allele_name\tnumber_of_found_in_database\tsamples_possessing_the_allele\n')
    for allele_name, list_sample in sorted(dict_report.items()):
        new_name = dict_old_name_2_new[allele_name]
        f_of.write(">" + new_name + '\n')
        f_of.write(dict_fasta[allele_name] + '\n')

        f_or.write(new_name + '\t' + str(len(list_sample)) + '\t' + ','.join(list_sample) + '\n')
    f_of.close()
    f_or.close()

