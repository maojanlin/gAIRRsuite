import argparse
import pickle
import os
import numpy as np
#from parse_contig_realign import mark_edit_region, variant_link_graph, haplotyping_link_graph, output_contig_correction
from parse_contig_realign import variant_link_graph, output_contig_correction, parse_CIGAR, parse_MD, trim_dict, find_double_pos, get_farthest_ext
from utils import get_reverse_complement
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'sam file of reads realign to contig'
    )
    parser.add_argument(
        '-fc', '--fn_cluster_contig',
        help = 'cropped contig file, corrected or not'
    )
    
    parser.add_argument(
        '-for', '--fo_report',
        help = 'output report file'
    )
    parser.add_argument(
        '-foc', '--fo_corrected_alleles',
        help = 'output corrected alleles fasta file'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



def cluster_separate(fn_cluster_contig, fn_sam):
    # dict_contig {}
    #  - keys: contig_name
    #  - values: [edit_histogram, cover_histogram, contig_SEQ, list_read_field[]]
    dict_contig = {}
    # dict_contig's initialization
    with open(fn_cluster_contig, 'r') as f_c:
        contig_name = ""
        contig_SEQ  = ""
        for line in f_c:
            if line[0] == '>':
                if contig_name != "":
                    dict_contig[contig_name] = [np.zeros(len(contig_SEQ) + 1), np.zeros(len(contig_SEQ) + 1), contig_SEQ, []]
                contig_name = line.strip()[1:].split()[0]
                contig_SEQ = ""
            else:
                contig_SEQ += line.strip()
        dict_contig[contig_name] = [np.zeros(len(contig_SEQ) + 1), np.zeros(len(contig_SEQ) + 1), contig_SEQ, []]

    with open(fn_sam, 'r') as f_r:
        read_name = ""
        read_SEQ = ""
        for line in f_r:
            if line[0] != '@':
                fields = line.split()
                if fields[2] == '*':
                    continue
                else:
                    contig_name = fields[2]
                    dict_contig[contig_name][3].append(fields)

    return dict_contig



def mark_edit_region(contig_name, contig_info, ignore_S=False):
    # contig_info = [edit_histogram, cov_histogram, contig_SEQ, list_read]
    edit_histogram = contig_info[0]
    cov_histogram  = contig_info[1]
    # list_read_info: [ (start_pos, end_pos, read_name, even_odd_flag, mis_region) ]
    list_read_info = []
    even_odd_flag = 1
    list_read_field = contig_info[3]
    for fields in list_read_field:
        read_name = fields[0]
        read_SEQ  = fields[9]
        cigar     = fields[5]
        sam_flag  = int(fields[1])
        # if the alignment is a supplementary alignment, pass, it does not matter the even odd
        # read BWA manual "Supplementary Alignment" for more information
        if sam_flag > 1024:
            continue

        S_flag = False
        number, operate = parse_CIGAR(cigar)
        if ignore_S and 'S' in cigar:
            if operate[0] == 'S':
                if number[0] >= len(read_SEQ)/15:
                    S_flag = True
            if operate[-1] == 'S':
                if number[-1] >= len(read_SEQ)/15:
                    S_flag = True
        # if cigar == '*', means alignment is bad, pass
        # if the read align to incorrect contigs, pass
        if cigar == '*' or contig_name != fields[2] or S_flag:
            # list_read_info.append((start_pos, end_pos, read_name, even_odd_flag, mis_region))
            list_read_info.append((0, 0, read_name, even_odd_flag, [], "", read_SEQ))
            if even_odd_flag == 1:
                even_odd_flag = 2
            else:
                even_odd_flag = 1
            continue

        edit_dist = int(fields[11].split(':')[2])  # NM:i:2 tag
        MD_tag    = fields[12].split(':')[2]       # MD:Z:38G2A20
        start_pos = int(fields[3])
        
        mis_region_MD = parse_MD(MD_tag)
        mis_region_MD = [ele + start_pos - 1 for ele in mis_region_MD] # change to ref coordinate

        mis_region_I = []   # insertion boundary region
        diff_len = 0        # len contribution of D, I, and S
        if 'I' in operate or 'D' in operate or 'S' in operate:
            idx_I = start_pos - 1 # index in reference
            for idx, op in enumerate(operate):
                if op == 'I':
                    diff_len -= number[idx]
                    mis_region_I.append(idx_I)
                    mis_region_I.append(idx_I+1)
                else:
                    if op == 'S':
                        diff_len -= number[idx]
                    else:
                        idx_I += number[idx]
                        if op == 'D':
                            diff_len += number[idx]
        
        end_pos = start_pos + len(fields[9]) + diff_len
        
        match_len = end_pos - start_pos
        mis_region_S = []
        recover_S_flag = False
        
        if operate[0] == 'S':
            left_S_len = min(number[0], start_pos-1)
            if left_S_len < match_len/10: # if S len is not too long, we accept it as mismatch
                mis_region_S = [pos for pos in range(start_pos-left_S_len,start_pos)]
                start_pos -= left_S_len
                operate[0] = 'M'
                if left_S_len != number[0]:
                    operate = ['S'] + operate 
                    number  = [number[0]-left_S_len] + number
                    number[1] = left_S_len
                recover_S_flag = True
        if operate[-1] == 'S':
            right_S_len = min(number[-1], len(cov_histogram)-end_pos)
            if right_S_len < match_len/10: # if S len is not to long, we accept it as mismatch
                mis_region_S += [pos for pos in range(end_pos,end_pos+right_S_len)]
                end_pos += right_S_len 
                operate[-1] = 'M'
                if right_S_len != number[-1]:
                    operate = operate + ['S']
                    number  = number  + [number[-1]-right_S_len]
                    number[-2] = right_S_len
                recover_S_flag = True
        if recover_S_flag:
            cigar = ""
            for cigar_id, element in enumerate(number):
                cigar += str(element)
                cigar += operate[cigar_id]
        
        #print(read_name + '\t', start_pos, end_pos)
        cov_histogram[start_pos:end_pos] += 1

        mis_region = mis_region_MD + mis_region_I + mis_region_S
        mis_region.sort()
        edit_histogram[mis_region] += 1
        
        # record the reads information
        list_read_info.append((start_pos, end_pos, read_name, even_odd_flag, mis_region, cigar, read_SEQ))
        if even_odd_flag == 1:
            even_odd_flag = 2
        else:
            even_odd_flag = 1

    return edit_histogram, cov_histogram, list_read_info


def haplotyping_link_graph(dict_link_graph, dict_var_weight, dict_link_outward, dict_link_inward, edit_region):
    # sort the potential variants on the interested site, can only use these variants bases
    list_pos_weight = []
    print("Trimming the significant bases at interested site:")
    print("Original site-base dict", dict_var_weight)
    for key in sorted(dict_var_weight.keys()):
        dict_part = dict_var_weight[key]
        trim_dict(dict_part, 10)
        list_pos_weight.append((key, sorted(dict_part.items(), key=lambda pair:pair[1], reverse=True)))

    print("Final site-base list:", list_pos_weight)
    eprint("#### max site-base variant #", max([len(ele[1]) for ele in list_pos_weight]))
    if list_pos_weight == []:
        print("There is no variant detected!")
        return [], []
    
    print("+++++++++++++++++++", "dict_link_graph", "+++++++++++++++++++")
    for key in sorted(dict_link_graph.keys()):
        print(key, dict_link_graph[key])
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    
    # initializing the haplotype list, the cursor, and the last_ext
    haplotype_0 = []        # record the (position, base) pair of the haplotype
    hap_cursor_0 = 0        # record the position got the linking information (still useless in this version)
    break_flag_0 = False    # the flag indicating of the haplotype is breaked
    haplotype_1 = []
    hap_cursor_1 = 0
    break_flag_1 = False
    pos_start_idx = 0
   
    # find the first variant site with two variants
    pos_start_idx, haplotype_0, haplotype_1, hap_cursor_0, hap_cursor_1 = find_double_pos(pos_start_idx, list_pos_weight, haplotype_0, haplotype_1, hap_cursor_0, hap_cursor_1)

    # haplotyping from list_pos_weight:
    for pos_idx in range(pos_start_idx, len(list_pos_weight)):
        pos_weight = list_pos_weight[pos_idx]
        position = pos_weight[0]
        list_pos_base = pos_weight[1]
        print("XXXXXXXXXXXXXX", position, "XXXXXXXXXXXXXXXX")

        # deal with haplotype_0's outward lin
        dict_outward_0 = {}
        if dict_link_outward.get(haplotype_0[hap_cursor_0]):
            dict_outward_0 = dict_link_outward[haplotype_0[hap_cursor_0]]
        trim_dict(dict_outward_0)
        if position > get_farthest_ext(dict_outward_0, haplotype_0[hap_cursor_0]):
            break_flag_0 = True
            eprint("Haplotype 0 has a break at", haplotype_0[hap_cursor_0], "to", position)
        print(dict_outward_0)
        # deal with haplotype_1's outward link
        print("--------------------")
        dict_outward_1 = {}
        if dict_link_outward.get(haplotype_1[hap_cursor_1]):
            dict_outward_1 = dict_link_outward[haplotype_1[hap_cursor_1]]
        trim_dict(dict_outward_1)
        if position > get_farthest_ext(dict_outward_1, haplotype_1[hap_cursor_1]):
            break_flag_1 = True
            eprint("Haplotype 1 has a break at", haplotype_1[hap_cursor_1], "to", position)
        print(dict_outward_1)
        # deal with position's inward link
        print("--------------------")
        dict_inward_0 = {}
        if dict_link_inward.get((position, list_pos_base[0][0])):
            dict_inward_0 = dict_link_inward[(position, list_pos_base[0][0])]
        trim_dict(dict_inward_0)
        print(dict_inward_0)
        #print(dict_link_graph[(position, list_pos_base[1][0])])
        if len(list_pos_base) > 1:
            print("--------------------")
            dict_inward_1 = {}
            if dict_link_inward.get((position, list_pos_base[1][0])):
                dict_inward_1 = dict_link_inward[(position, list_pos_base[1][0])]
            trim_dict(dict_inward_1)
            print(dict_inward_1)

        connect_info_0 = None
        connect_info_1 = None
        # There must be at least one kind of base in the position
        for (outward_key, weight) in sorted(dict_outward_0.items(), key=lambda pair:pair[1], reverse=True):
            if dict_inward_0.get(outward_key):
                print("Potential Connect: ", outward_key, 0, 0)
                connect_info_0 = (dict_outward_0[outward_key], (position, outward_key[1][1]))
                break
        for (outward_key, weight) in sorted(dict_outward_1.items(), key=lambda pair:pair[1], reverse=True):
            if dict_inward_0.get(outward_key):
                print("Potential Connect: ", outward_key, 1, 0)
                connect_info_1 = (dict_outward_1[outward_key], (position, outward_key[1][1]))
                break
        
        # if there are two variants in the position
        if len(list_pos_base) > 1:
            for (outward_key, weight) in sorted(dict_outward_0.items(), key=lambda pair:pair[1], reverse=True):
                if dict_inward_1.get(outward_key):
                    print("Potential Connect: ", outward_key, 0, 1)
                    if connect_info_0 == None or connect_info_0[0] < weight:
                        connect_info_0 = (dict_outward_0[outward_key], (position, outward_key[1][1]))
                        break
            for (outward_key, weight) in sorted(dict_outward_1.items(), key=lambda pair:pair[1], reverse=True):
                if dict_inward_1.get(outward_key):
                    print("Potential Connect: ", outward_key, 1, 1)
                    if connect_info_1 == None or connect_info_1[0] < weight:
                        connect_info_1 = (dict_outward_1[outward_key], (position, outward_key[1][1]))
                        break
        
        # the case that two haplotypes may collapse into one
        if connect_info_0 and connect_info_1:
            if connect_info_0[1] == connect_info_1[1]: # two haplotypes are collapsed
                record_info_0 = [connect_info_0[1]]
                record_info_1 = [connect_info_1[1]]
                for redouble_idx in range(pos_idx, len(list_pos_weight)):
                    rd_pos_weight = list_pos_weight[redouble_idx]
                    rd_position = rd_pos_weight[0]
                    rd_list_pos_base = rd_pos_weight[1]
                    if(len(rd_list_pos_base)) >= 2: # if there are two variants at the site
                        # call the potential connections
                        last_info_0 = haplotype_0[hap_cursor_0]
                        last_info_1 = haplotype_1[hap_cursor_1]
                        dict_info_0 = dict_link_graph[last_info_0]
                        dict_info_1 = dict_link_graph[last_info_1]
                        # connect them
                        rd_info_0 = None
                        rd_info_1 = None
                        for rd_link_info, rd_weight in sorted(dict_info_0.items(), key=lambda pair:pair[1], reverse=True):
                            variant_flag = False
                            for info_pair in rd_link_info[1]:
                                tmp_rd_info = []
                                if info_pair == connect_info_0[1]:
                                    variant_flag = True
                                    tmp_rd_info.append(info_pair)
                                if variant_flag:
                                    tmp_rd_info.append(info_pair)
                                    if info_pair[0] == rd_position:
                                        rd_info_0 = tmp_rd_info
                                        break
                            if rd_info_0:
                                break
                            
                        for rd_link_info, rd_weight in sorted(dict_info_1.items(), key=lambda pair:pair[1], reverse=True):
                            for info_pair in rd_link_info[1]:
                                tmp_rd_info = []
                                if info_pair == connect_info_1[1]:
                                    variant_flag = True
                                    tmp_rd_info.append(info_pair)
                                if variant_flag:
                                    tmp_rd_info.append(info_pair)
                                    if info_pair[0] == rd_position:
                                        rd_info_1 = tmp_rd_info
                                        break
                            if rd_info_1:
                                break

                        print("connect_info_0", record_info_0)
                        print("connect_info_1", record_info_1)
                        print("rd_info_0", rd_info_0)
                        print("rd_info_1", rd_info_1)
                        if rd_info_0:
                            record_info_0 += rd_info_0
                        if rd_info_1:
                            record_info_1 += rd_info_1
                        if rd_info_0 != rd_info_1:
                            if rd_info_0:
                                pass
                            else:
                                break_flag_0 = True
                            if rd_info_1:
                                pass
                            else:
                                break_flag_1 = True
                            break
                haplotype_0 += record_info_0
                hap_cursor_0 += len(record_info_0)
                haplotype_1 += record_info_1
                hap_cursor_1 += len(record_info_1)
                print("Crossing the single base variant site...")
                continue



        # update the nodes if the connection is found
        if connect_info_0:    
            haplotype_0.append(connect_info_0[1])
            hap_cursor_0 += 1
            if break_flag_1 and len(list_pos_base) >1:
                for idx in range(2):
                    potential_base = list_pos_base[idx][0] 
                    if potential_base != connect_info_0[1][1]:
                        eprint("Link rebuilt on Haplotype 1 at", haplotype_1[hap_cursor_1] , "to", position)
                        haplotype_1.append((position, potential_base))
                        hap_cursor_1 += 1
                        break_flag_1 = False
                        break
        if connect_info_1:    
            haplotype_1.append(connect_info_1[1])
            hap_cursor_1 += 1
            if break_flag_0 and len(list_pos_base) >1:
                for idx in range(2):
                    potential_base = list_pos_base[idx][0] 
                    if potential_base != connect_info_1[1][1]:
                        eprint("Link rebuilt on Haplotype 0 at", haplotype_0[hap_cursor_0] , "to", position)
                        haplotype_0.append((position, potential_base))
                        hap_cursor_0 += 1
                        break_flag_0 = False
                        break

        if break_flag_0 and break_flag_1:
            eprint("BREAKING LINKS FOR BOTH HAPLOTYPE AT", position, "!!!!")
            eprint("Breaking links cannot resloved, we guess...")
            haplotype_0.append((position, list_pos_base[0][0]))
            hap_cursor_0 += 1
            if len(list_pos_base) > 1:
                haplotype_1.append((position, list_pos_base[1][0]))
                hap_cursor_1 += 1
    
    print(haplotype_0)
    print(haplotype_1)
    return haplotype_0, haplotype_1





if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    fn_cluster_contig = args.fn_cluster_contig
    
    fo_report = args.fo_report
    fo_corrected_alleles = args.fo_corrected_alleles
    
    dict_contig = cluster_separate(fn_cluster_contig, fn_sam)

    for contig_name, contig_info in sorted(dict_contig.items()):
        #parse the sam file and generate
        edit_histogram, cov_histogram, list_read_info = mark_edit_region(contig_name, contig_info)
        
        #determine the region contains alternative flanking region
        edit_region = []
        for idx, ele in enumerate(edit_histogram):
            print(str(idx) + ':\t' + str(cov_histogram[idx])  + '\t' + str(ele))
            if ele > cov_histogram[idx]/4:
                edit_region.append(idx)

        print(contig_name, edit_region)
        
        contig_SEQ = dict_contig[contig_name][2]
        interest_region = "0-" + str(len(contig_SEQ))
        interest_edit_region = edit_region
        if interest_edit_region != [] and min(cov_histogram[1:]) > 20:
            print("=========== allele correction ==============")
            eprint("CORRECT", contig_name.split('|')[1], min(cov_histogram[1:]), interest_edit_region)
            dict_link_graph, dict_var_weight, dict_link_outward, dict_link_inward = variant_link_graph(interest_edit_region, list_read_info)
            haplotype_0, haplotype_1 = haplotyping_link_graph(dict_link_graph, dict_var_weight, dict_link_outward, dict_link_inward, interest_region)
            #output_contig_correction(contig_SEQ, region_st, region_ed, haplotype_0, haplotype_1, contig_name, corrected_contig_output_file)
            diff_0 = sum([max(len(ele[1])-1,0) for ele in haplotype_0 ])
            diff_1 = sum([max(len(ele[1])-1,0) for ele in haplotype_1 ])
            len_insert = max(diff_0, diff_1)

            output_contig_correction(contig_SEQ, 0, len(contig_SEQ)+len_insert, haplotype_0, haplotype_1, contig_name, fo_corrected_alleles, "/novel")
        elif interest_edit_region != []:
            eprint("DDeficient", contig_name.split('|')[1], min(cov_histogram[1:]), interest_edit_region)
            print("=== cov not efficient:", min(cov_histogram[1:]), "=======")
        else:
            eprint("No variant", contig_name.split('|')[1], min(cov_histogram[1:]), interest_edit_region)
            print("============ No novel allele ===============")
     

