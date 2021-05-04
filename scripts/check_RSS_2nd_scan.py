import argparse
from utils import fasta_to_dict, parse_CIGAR
from check_RSS_1st_scan import get_ref_len

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'input sam file that hepNona_file align to flanking haplotypes'
    )
    parser.add_argument(
        '-ff', '--fn_flanking',
        help = 'input fasta file that contains flanking haplotypes'
    )
    parser.add_argument(
        '-gtype', '--gene_type',
        help = '\"TCRV\", \"BCRV\" for V alleles and \"TCRJ\" for J alleles'
    )
    parser.add_argument(
        '-fo', '--fo_output',
        help = 'output csv file'
    )
    parser.add_argument(
        '-fos', '--fo_summary',
        help = 'output summary file'
    )
    parser.add_argument(
        '-foh', '--fo_database',
        help = 'output database csv file'
    )
    parser.add_argument(
        '-fon', '--fo_novel_fasta',
        help = 'output novel RSS fasta file'
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
    CIGAR = fields[5]
    list_num, list_operand = parse_CIGAR(CIGAR)
    ref_len = get_ref_len(list_num, list_operand)
    if len(fields) >= 11:
        NM = int(fields[11].split(':')[2])
    return (strend, position, NM, ref_len)


def parse_sam_with_lenInfo(fn_sam):
    # dict_haplotype_hepNona = {}
    # -  key: haplotype_name
    # -  values: [haplotype_len, set(heptamer info), set(nonamer info)]
    #        - nonamer info: (strend, position, ref_len)
    dict_haplotype_hepNona = {}

    f_n = open(fn_sam)
    for line in f_n:
        fields = line.split()
        if fields[0] == "@SQ":
            haplotype_name = fields[1][3:]
            haplotype_len  = int(fields[2][3:])
            dict_haplotype_hepNona[haplotype_name] = [haplotype_len, set(), set()]
        elif ('HEPTAMER' in fields[0]):
            haplotype_name = fields[2]
            if haplotype_name == '*':
                continue
            info_tuple = hep_nano_info(fields)
            if info_tuple[0] and (info_tuple[2] <= 2):
                dict_haplotype_hepNona[haplotype_name][1].add( (info_tuple[0], info_tuple[1], info_tuple[3]) )
        elif ('NONAMER' in fields[0]):
            haplotype_name = fields[2]
            if haplotype_name == '*':
                continue
            info_tuple = hep_nano_info(fields)
            if info_tuple[0] and (info_tuple[2] <= 2):
                dict_haplotype_hepNona[haplotype_name][2].add( (info_tuple[0], info_tuple[1], info_tuple[3]) )
    f_n.close()
    return dict_haplotype_hepNona


def select_best_RSS(interested_list, gene_type):
    """ranks are list as following: perfect, shift perfect, +/-1, shift +/-1"""
    ranks = [[],[],[],[]]
    perfect_len = 23
    if gene_type == 1:
        pass
    elif gene_type == 0:
        perfect_len = -12
    else:
        print("WANRING: not support IG yet...")
    
    for pair in interested_list:
        if (pair[0] == 0) and (pair[1] == perfect_len):
            ranks[0].append(pair)
        elif (-10 <= pair[0] <= 10) and (pair[1] == perfect_len):
            ranks[1].append(pair)
        elif pair[0] == 0:
            ranks[2].append(pair)
        else:
            ranks[3].append(pair)
    for idx, sub_rank in enumerate(ranks):
        if len(sub_rank) > 0:
            return idx, sub_rank
    return 4, [(-1,-1,-1,-1,"")]



if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    fn_flanking = args.fn_flanking
    fo_output = args.fo_output
    fo_summary = args.fo_summary
    fo_database = args.fo_database
    fo_novel_fasta = args.fo_novel_fasta
    gene_type = args.gene_type
    if "TCR" or "tcr" in gene_type:
        if 'V' in gene_type or 'v' in gene_type:
            gene_type = 1
        else:
            gene_type = 0
    elif "BCR" or "bcr" in gene_type:
        if 'V' in gene_type or 'v' in gene_type:
            gene_type = 3
        else:
            gene_type = 2
    else:
        print("WARNING: ERROR in gene_type!!")

    dict_flank = fasta_to_dict(fn_flanking)
    dict_haplotype_hepNona = parse_sam_with_lenInfo(fn_sam)

    # dict_list_fault:
    # - key: gene_name
    # - value: set( print_pair_info )
    list_fault = [{},{},{},{},{}]
    list_recover_hap = []
    list_missing_hap = []
    set_flank  = set()
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
                #if ( 11 <= abs(spacer_len) <= 13 ) or ( 22 <= abs(spacer_len) <= 24): # make sure within 12/23
                if gene_type == 1:
                    if (20 <= spacer_len <= 30):
                        normalized_pos = position_hep + 199 - haplotype_len
                        haplotype_SEQ = dict_flank[haplotype_name]
                        RSS_len = spacer_len + hep_info_tuple[2] + nona_info_tuple[2]
                        RSS_SEQ = haplotype_SEQ[position_hep-1:position_hep-1+RSS_len]
                        if normalized_pos > -100:
                            interested_pair.append( (normalized_pos, spacer_len, hep_info_tuple[0], nona_info_tuple[0], RSS_SEQ) )
                elif gene_type == 0:
                    if (-13 <= spacer_len <= -7):
                        normalized_pos = position_hep - 194
                        haplotype_SEQ = dict_flank[haplotype_name]
                        RSS_len = abs(spacer_len) + hep_info_tuple[2] + nona_info_tuple[2]
                        RSS_SEQ = haplotype_SEQ[position_nona-1:position_nona-1+RSS_len]
                        if normalized_pos < 60:
                            interested_pair.append( (position_hep - 194, spacer_len, hep_info_tuple[0], nona_info_tuple[0], RSS_SEQ) )
                else:
                    print("WARNING: This program do not support IG yet...")
        rank, best_RSS = select_best_RSS(interested_pair, gene_type)
        ext_seq = dict_flank[haplotype_name]
        
        # setup for allele level group analysis
        allele_name = haplotype_name[:haplotype_name.rfind('/')]
        f_o.write(haplotype_name + ',' + str(haplotype_len) + ',' + str(rank) + ',')
        print_pair_info = ""
        for pair in best_RSS:
            print_pair_info += ("(pos: " + str(pair[0]) + ', spacer: ' + str(pair[1]) + ', SEQ: ' + str(pair[4]) + "),")
            keep_info = (pair[0], pair[1], pair[4])
            if list_fault[rank].get(allele_name):
                if keep_info in list_fault[rank][allele_name]:
                    pass
                else:
                    list_fault[rank][allele_name].add(keep_info)
            else:
                list_fault[rank][allele_name] = {keep_info}
        f_o.write(print_pair_info + "\n")
        
        # setup for database data
        if rank == 4:
            list_missing_hap.append([haplotype_name,allele_functionality(allele_name)]) # pair[4] is the SEQ
        else:
            for pair in best_RSS:
                list_recover_hap.append([haplotype_name,allele_functionality(allele_name), pair[4]])
    f_o.close()

    f_h = open(fo_database, 'w')
    f_h.write("Allele_without_RSS:" + str(len(list_missing_hap)) + '\n')
    f_h.write("flanking_sequence_name,gene_functionality,RSS_id,sequence\n")
    for haplotype_name, functionality in sorted(list_missing_hap):
        f_h.write(haplotype_name + ',' + functionality + '\n')
    f_h.write("Allele_with_novel_RSS:" + str(len(list_recover_hap)) + '\n')
    f_h.write("flanking_sequence_name,gene_functionality,RSS_id,sequence\n")
    dict_allele_count = {}
    for pair_info in sorted(list_recover_hap):
        haplotype_name = pair_info[0]
        allele_name = haplotype_name[:haplotype_name.rfind('/')]
        functionality  = pair_info[1]
        SEQ            = pair_info[2] 
        if dict_allele_count.get(allele_name):
            if SEQ in dict_allele_count[allele_name]:
                pass
            else:
                dict_allele_count[allele_name].add(SEQ)
            count_id = len(dict_allele_count[allele_name])
            f_h.write(haplotype_name + ',' + functionality + ',' + allele_name + '_RSS*n' + str(count_id).zfill(2) + ',' + SEQ + '\n')
        else:
            dict_allele_count[allele_name] = {SEQ}
            count_id = 1
            f_h.write(haplotype_name + ',' + functionality + ',' + allele_name + '_RSS*n' + str(count_id).zfill(2) + ',' + SEQ + '\n')
    
    f_s = open(fo_summary, 'w')
    rank_property = ["perfect", "shift (10) perfect", "+/-1", "shift (10) +/-1", "Not found"]
    # ========== print alleles with missing or novel RSS ==========
    f_s.write("Heptamer/Nonamer cannot be found: " + str(len(list_fault[4])) + '\n')
    for allele_name, print_pair_info in sorted(list_fault[4].items()):
        f_s.write("\t" + allele_name + '\t' + allele_functionality(allele_name) + '\t' + str(print_pair_info) + '\n')
        #f_h.write(allele_name + ',' + allele_functionality(allele_name) + '\n')
    # ========== print alleles with +/- 1 spacer ==========
    set_violate = set(list_fault[2].keys()) | set(list_fault[3].keys())
    set_normal  = set(list_fault[0].keys()) | set(list_fault[1].keys())
    f_s.write("Alleles with abnormal spacer: " + str(len(set_violate)) + '\n')
    for allele_name in sorted(set_violate):
        if list_fault[2].get(allele_name):
            print_pair_info = list_fault[2][allele_name]
        else:
            print_pair_info = list_fault[3][allele_name]
        f_s.write("\t" + allele_name + '\t' + allele_functionality(allele_name) + '\t' + str(print_pair_info) + '\n')
        #for pair_info in print_pair_info:
        #    f_h.write(allele_name + ',' + allele_functionality(allele_name) + ',,' + pair_info[2] + '\n')

    # ========== print alleles with normal or shift RSS ==========
    f_s.write("Alleles with 12/23 spacer: " + str(len(set_normal)) + '\n')
    for allele_name in sorted(set_normal):
        if list_fault[0].get(allele_name):
            print_pair_info = list_fault[0][allele_name]
        else:
            print_pair_info = list_fault[1][allele_name]
        f_s.write("\t" + allele_name + '\t' + allele_functionality(allele_name) + '\t' + str(print_pair_info) + '\n')
        #for pair_info in print_pair_info:
        #    f_h.write(allele_name + ',' + allele_functionality(allele_name) + ',,' + pair_info[2] + '\n')
    f_s.close()
    f_h.close()

    f_n = open(fo_novel_fasta, 'w')
    for rank in (0,1,2,3):
        for allele_name, set_info in sorted(list_fault[rank].items()):
            for idx, info in enumerate(sorted(set_info)):
                f_n.write('>' + allele_name + '/RSS*n' + str(idx+1).zfill(2) + '\n')
                f_n.write(info[2] + '\n')
    f_o.close()
