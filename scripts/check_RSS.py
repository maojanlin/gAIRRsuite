import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'input sam file that hepNona_file align to flanking haplotypes'
    )
    parser.add_argument(
        '-gtype', '--gene_type',
        help = '\"TCRV\" for V alleles and \"TCRJ\" for J alleles'
    )
    parser.add_argument(
        '-fo', '--fo_output',
        help = 'output csv file'
    )
    parser.add_argument(
        '-fos', '--fo_summary',
        help = 'output summary file'
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



def hep_nano_info(fields):
    strend = True
    if int(fields[1]) % 32 >= 16:
        strend = False
    position = int(fields[3])
    if len(fields) >= 11:
        NM = fields[11]
    return (strend, position, NM)


def parse_sam_with_lenInfo(fn_sam):
    # dict_haplotype_hepNona = {}
    # -  key: haplotype_name
    # -  values: [haplotype_len, set(heptamer info), set(nonamer info)]
    #        - nonamer info: (strend, position, NM:i)
    dict_haplotype_hepNona = {}

    f_n = open(fn_sam)
    for line in f_n:
        fields = line.split()
        if fields[0] == "@SQ":
            haplotype_name = fields[1][3:]
            haplotype_len  = int(fields[2][3:])
            dict_haplotype_hepNona[haplotype_name] = [haplotype_len, set(), set()]
        elif ('HEPTAMER' in fields[0]):
            heplotype_name = fields[2]
            info_tuple = hep_nano_info(fields)
            if info_tuple[0] and (info_tuple[2] == "NM:i:0"):
                dict_haplotype_hepNona[heplotype_name][1].add( info_tuple )
        elif ('NONAMER' in fields[0]):
            heplotype_name = fields[2]
            info_tuple = hep_nano_info(fields)
            if info_tuple[0] and (info_tuple[2] == "NM:i:0"):
                dict_haplotype_hepNona[heplotype_name][2].add( info_tuple )
    f_n.close()
    return dict_haplotype_hepNona


def select_best_RSS(interested_list):
    ranks = [[],[],[],[]]
    for pair in interested_list:
        if (pair[0] == 0) and (abs(pair[1]) == 12 or abs(pair[1]) == 23):
            ranks[0].append(pair)
        elif pair[0] == 0:
            ranks[1].append(pair)
        elif (-10 <= pair[0] <= 10) and (pair[1] == 12 or pair[1] == 23):
            ranks[2].append(pair)
        else:
            ranks[3].append(pair)
    for idx, sub_rank in enumerate(ranks):
        if len(sub_rank) > 0:
            return idx, sub_rank
    return 4, interested_list



if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    fo_output = args.fo_output
    fo_summary = args.fo_summary
    gene_type = args.gene_type
    if 'V' in gene_type or 'v' in gene_type:
        gene_type = 1
    else:
        gene_type = 0

    dict_haplotype_hepNona = parse_sam_with_lenInfo(fn_sam)

    list_fault = []
    list_bad   = []
    list_num   = [0,0,0,0,0]
    f_o = open(fo_output, 'w')
    f_o.write("Extended_allele_name,len,rank,RSS_pairs\n")
    for haplotype_name, info_pair in sorted(dict_haplotype_hepNona.items()):
        haplotype_len = info_pair[0]
        heptamer_info = sorted(info_pair[1])
        nonamer_info  = sorted(info_pair[2])
        interested_pair = []
        for hep_info_tuple in heptamer_info:
            position_hep = hep_info_tuple[1]
            for nona_info_tuple in nonamer_info:
                position_nona = nona_info_tuple[1]
                spacer_len = position_nona - position_hep
                # manage the distance between nearest ends of heptamer and nonamer
                if spacer_len > 7:
                    spacer_len -= 7
                elif spacer_len < 9:
                    spacer_len += 9
                else:
                    break
                if ( 11 <= abs(spacer_len) <= 13 ) or ( 22 <= abs(spacer_len) <= 24): # make sure within 12/23
                    if gene_type == 1:
                        interested_pair.append( (position_hep + 199 - haplotype_len, spacer_len, hep_info_tuple[0], nona_info_tuple[0]) )
                    else:
                        interested_pair.append( (position_hep - 194, spacer_len, hep_info_tuple[0], nona_info_tuple[0]) )
        rank, best_RSS = select_best_RSS(interested_pair)
        
        f_o.write(haplotype_name + ',' + str(haplotype_len) + ',' + str(rank) + ',')
        print_pair_info = ""
        for pair in best_RSS:
            print_pair_info += ("(" + str(pair[0]) + ' ' + str(pair[1]) + "),")
        f_o.write(print_pair_info + "\n")
        
        list_num[rank] += 1
        if 1 <= rank <= 3:
            list_fault.append((haplotype_name, print_pair_info))
        elif rank >= 4:
            list_bad.append((haplotype_name, print_pair_info))
    f_o.close()

    f_s = open(fo_summary, 'w')
    f_s.write("Extended_allele with perfect RSS: " + str(list_num[0]) + '\n')
    f_s.write("Extended_allele with fault RSS (+/-10 position or +/- 12/23 spacer): " + str(sum(list_num[1:4])) + '\n')
    for haplotype_name, print_pair_info in list_fault:
        f_s.write('\t' + haplotype_name + '\t' + print_pair_info + '\n')
    f_s.write("Extended_allele without known RSS: " + str(list_num[4]) + '\n')
    for haplotype_name, print_pair_info in list_bad:
        f_s.write('\t' + haplotype_name + '\t' + allele_functionality(haplotype_name) + '\t' + print_pair_info + '\n')
    f_s.close()
