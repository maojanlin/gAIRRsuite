import argparse
from utils import fasta_to_dict, parse_CIGAR, get_reverse_complement

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'input sam file that RSS_file align to flanking haplotypes'
    )
    parser.add_argument(
        '-ff', '--fn_flanking',
        help = 'input initial fasta file that contains flanking haplotypes'
    )
    parser.add_argument(
        '-fo', '--fo_output',
        help = 'output csv file'
    )
    parser.add_argument(
        '-fod', '--fo_detail',
        help = 'output detailed haplotype report file'
    )
    parser.add_argument(
        '-fos', '--fo_summary',
        help = 'output summary allele report file'
    )
    parser.add_argument(
        '-fonf', '--fo_novel_fasta',
        help = 'output novel RSS fasta file'
    )
    #parser.add_argument(
    #    '-fonc', '--fo_novel_csv',
    #    help = 'output novel RSS csv file'
    #)
    #parser.add_argument(
    #    '-fooc', '--fo_original_csv',
    #    help = 'output original RSS csv file'
    #)
    parser.add_argument(
        '-fomf', '--fo_missing_fasta',
        help = 'output missing RSS fasta file for next stage analysis'
    )
    args = parser.parse_args()
    return args


set_F_ORF = {'TRBJ2-7', 'TRBV7-3'}
set_F_P   = {"TRAJ8", "TRAV29/DV5", "TRAV3", "TRAV35", "TRBV10-1", "TRBV16", "TRBV30", "TRBV7-4"}
set_ORF = {'TRAJ1', 'TRAJ19', 'TRAJ2' , 'TRAJ25', 'TRAJ58', 'TRAJ59', 'TRAJ61', 'TRBJ2-2P', 'TRBV17', 'TRBV20/OR9-2', 'TRBV23-1', 'TRBV23/OR9-2', 'TRBV24/OR9-2', 'TRBV29/OR9-2', 'TRBV5-3', 'TRBV5-7', 'TRBV6-7', 'TRBV7-1', 'TRGV1', 'TRGV10', 'TRGV11'}
set_P = {'TRAJ51', 'TRAJ55', 'TRAJ60', 'TRAV11', 'TRAV11-1', 'TRAV14-1', 'TRAV15', 'TRAV28', 'TRAV31', 'TRAV32', 'TRAV33', 'TRAV37', 'TRAV46', 'TRAV8-5', 'TRAV8-6-1', 'TRAV8-7', 'TRAVA', 'TRAVB', 'TRAVC', 'TRBV1', 'TRBV12-1', 'TRBV12-2', 'TRBV21-1', 'TRBV21/OR9-2', 'TRBV22-1', 'TRBV22/OR9-2', 'TRBV23/OR9-2', 'TRBV24/OR9-2', 'TRBV25/OR9-2', 'TRBV26', 'TRBV26/OR9-2', 'TRBV3-2', 'TRBV5-2', 'TRBV7-5', 'TRBV8-1', 'TRBV8-2', 'TRBVA', 'TRBVA/OR9-2', 'TRBVB', 'TRBVC', 'TRGV5P', 'TRGV6', 'TRGV7', 'TRGVA', 'TRGVB'}


def allele_functionality(allele_name):
    gene_name = allele_name.split('*')[0]
    if gene_name in set_ORF:
        return "ORF"
    elif gene_name in set_P:
        return "Pseudo"
    elif gene_name in set_F_P:
        return "F P"
    elif gene_name in set_F_ORF:
        return "F_ORF"
    else:
        return "Func"


def get_ref_len(CIGAR_num, CIGAR_operand):
    ref_len = 0
    for idx, operand in enumerate(CIGAR_operand):
        if (operand == 'M') or (operand == 'H') or (operand == 'S') or (operand == 'I'):
            ref_len += CIGAR_num[idx]
        elif operand == 'D':
            pass
        else:
            print("WARNING: CIGAR operand error!!!")
            return None
    return ref_len


def parse_sam(fn_sam, dict_ref):
    # dict_haplotype_hepNona = {}
    # -  key: haplotype_name
    # -  values: [haplotype_len, set(heptamer info), set(nonamer info)]
    #        - nonamer info: (strend, position, NM:i)
    dict_haplotype_RSS = {}

    f_n = open(fn_sam)
    for line in f_n:
        fields = line.split()
        if fields[0] == "@SQ":
            continue
            haplotype_name = fields[1][3:]
            allele_name = haplotype_name[:haplotype_name.rfind('/')]
            haplotype_len  = int(fields[2][3:])
            dict_haplotype_RSS[haplotype_name] = [haplotype_len, set(), set()]
        elif ('RSS' in fields[0]):
            RSS_name = fields[0]
            haplotype_name = fields[2]
            if haplotype_name == '*':
                continue
            strend = True
            if int(fields[1]) % 32 >= 16:
                strend = False
            position = int(fields[3])
            NM = int(fields[11].split(':')[2])
            (list_num, list_operand) = parse_CIGAR(fields[5])
            if (list_operand[0] == 'H') or (list_operand[0] == 'S'):
                position -= list_num[0]
                NM += list_num[0]
            if (list_operand[-1] == 'H') or (list_operand[-1] == 'S'):
                NM += list_num[-1]

            ref_len = get_ref_len(list_num, list_operand)
            RSS_SEQ = dict_ref[haplotype_name][position-1:position+ref_len-1]

            if dict_haplotype_RSS.get(haplotype_name):
                dict_haplotype_RSS[haplotype_name].append((position, NM, RSS_name, RSS_SEQ))
            else:
                dict_haplotype_RSS[haplotype_name] = [(position, NM, RSS_name, RSS_SEQ)]
    
    f_n.close()
    return dict_haplotype_RSS


def merge_haplotype_to_allele(list_RSS):
    dict_allele_RSS = {}
    for info_pair in list_RSS:
        haplotype_name = info_pair[0]
        best_matched_RSS = info_pair[1]
        allele_name = haplotype_name[:haplotype_name.rfind('/')]
        if dict_allele_RSS.get(allele_name):
            another_best_match_RSS = dict_allele_RSS[allele_name]
            if best_matched_RSS == another_best_match_RSS:
                continue
            else:
                for RSS_info_pair in best_matched_RSS:
                    exist_flag = False
                    RSS_SEQ = RSS_info_pair[3]
                    for another_RSS_info_pair in another_best_match_RSS:
                        if RSS_SEQ.lower() == another_RSS_info_pair[3].lower():
                            exist_flag = True
                            break
                    if exist_flag == False:
                        dict_allele_RSS[allele_name].append(RSS_info_pair)
        else:
            dict_allele_RSS[allele_name] = best_matched_RSS
    return dict_allele_RSS






if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    fn_flanking = args.fn_flanking
    fo_detail   = args.fo_detail
    fo_summary  = args.fo_summary
    fo_output   = args.fo_output
    fo_novel_fasta   = args.fo_novel_fasta
    #fo_novel_csv     = args.fo_novel_csv
    #fo_original_csv  = args.fo_original_csv
    fo_missing_fasta = args.fo_missing_fasta

    dict_flank = fasta_to_dict(fn_flanking)
    set_allele_name = {v[:v.rfind('/')] for v in sorted(dict_flank.keys())}
    dict_haplotype_RSS = parse_sam(fn_sam, dict_flank)

    list_perfect_RSS = []
    list_novel_RSS = []
    list_no_RSS = []
    #for haplotype_name in sorted(set_allele_name):
    for haplotype_name in sorted(dict_flank.keys()):
        if dict_haplotype_RSS.get(haplotype_name):
            matched_RSS = dict_haplotype_RSS[haplotype_name]
            best_matched_RSS = []
            # match_info: (position, NM, RSS_name)
            for match_info in matched_RSS:
                collision_flag = False
                for idx, best_info in enumerate(best_matched_RSS):
                    if abs(match_info[0] - best_info[0]) < 15:
                        collision_flag = True
                        if match_info[1] < best_info[1]:
                            best_matched_RSS[idx] = match_info
                if collision_flag == False:
                    best_matched_RSS.append(match_info)
            
            if best_matched_RSS[0][1] == 0:
                list_perfect_RSS.append((haplotype_name, best_matched_RSS))
            else:
                list_novel_RSS.append((haplotype_name, best_matched_RSS))
        else:
            list_no_RSS.append((haplotype_name))

    # ========== output missing RSS's flanking sequences for next stage analysis ===========
    f_om = open(fo_missing_fasta, 'w')
    for info_pair in list_no_RSS:
        haplotype_name = info_pair
        haplotype_SEQ = dict_flank[haplotype_name]
        f_om.write('>' + haplotype_name + '\n')
        f_om.write(haplotype_SEQ)
        f_om.write('\n')
    f_om.close()

    # ========== output the haplotype detailed file ===============
    f_o = open(fo_detail, 'w')
    f_o.write("RSS not found: " + str(len(list_no_RSS)) + '\n')
    for element in list_no_RSS:
        f_o.write('\t' + element + '\t' + allele_functionality(element) + '\n')
    f_o.write("Alleles with novel RSS: " + str(len(list_novel_RSS)) + '\n')
    for element in list_novel_RSS:
        f_o.write('\t' + str(element) + '\t' + allele_functionality(element[0])  + '\n')
    f_o.write("Alleles with perfect RSS: " + str(len(list_perfect_RSS)) + '\n')
    for element in list_perfect_RSS:
        f_o.write('\t' + str(element) + '\t' + allele_functionality(element[0])  + '\n')
    f_o.close()

    # ============ output the novel RSS sequences ==================
    dict_allele_novel_RSS = merge_haplotype_to_allele(list_novel_RSS)
    f_of = open(fo_novel_fasta, 'w')
    for allele_name, list_RSS_info in sorted(dict_allele_novel_RSS.items()):
        for idx, RSS_info in enumerate(list_RSS_info):
            f_of.write('>' + allele_name + '/RSS*n' + str(idx).zfill(2) + '\n')
            f_of.write(RSS_info[3] + '\n')
    f_of.close()

    # ============ output all RSS based on alleles into summary file ==========
    set_allele_missing_RSS = {name[:name.rfind('/')] for name in list_no_RSS}
    dict_allele_perfect_RSS = merge_haplotype_to_allele(list_perfect_RSS)
    f_o = open(fo_summary, 'w')
    f_o.write("RSS not found: " + str(len(set_allele_missing_RSS)) + '\n')
    for allele_name in sorted(set_allele_missing_RSS):
        f_o.write('\t' + allele_name + '\t' + allele_functionality(allele_name) + '\n')
    f_o.write("Alleles with novel RSS: " + str(len(dict_allele_novel_RSS)) + '\n')
    for allele_name, info in sorted(dict_allele_novel_RSS.items()):
        f_o.write('\t' + allele_name + '\t' + allele_functionality(allele_name) + '\t' + str(info) + '\n')
    f_o.write("Alleles with perfect RSS: " + str(len(dict_allele_perfect_RSS)) + '\n')
    for allele_name, info in sorted(dict_allele_perfect_RSS.items()):
        f_o.write('\t' + allele_name + '\t' + allele_functionality(allele_name) + '\t' + str(info) + '\n')
    f_o.close()


