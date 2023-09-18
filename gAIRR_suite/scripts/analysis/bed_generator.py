import pysam
import argparse

TRA_chain = {'TRA', 'TRD'}
TRB_chain = {'TRB'}
TRG_chain = {'TRG'}
IGH_chain = {'IGH'}
IGL_chain = {'IGL'}
IGK_chain = {'IGK'}



def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-rl', '--list_report', nargs='+', default=[], help='list of the report files of all the TR or IG')
    parser.add_argument('-o', '--output', help='the output report file')
    args = parser.parse_args(arguments)
    
    list_report = args.list_report
    fn_out      = args.output

    # read all the data files
    list_data = [[],[]]
    for report_name in list_report:
        idx = 0
        f = open(report_name, 'r')
        for line in f:
            fields = line.strip().split(',')
            if len(fields) < 2:
                idx = (idx + 1)%2
                continue
            elif len(fields) > 5:
                allele_name = fields[0] + '(' + fields[4][3:] + ',h' + fields[5][4:] + ')'
            elif len(fields) > 4:
                allele_name = fields[0] + '(' + fields[4][3:] + ')'
            else:
                allele_name = fields[0]
            list_data[idx].append((fields[1],int(fields[2]),int(fields[3]),allele_name))

    # sort the data files
    dict_chain = [{},{}]
    dict_chain[0] = {'TR/OR':[], 'IG/OR':[], 'TRA/D':[], 'TRB':[], 'TRG':[], 'IGH':[], 'IGL':[], 'IGK':[]}
    dict_chain[1] = {'TR/OR':[], 'IG/OR':[], 'TRA/D':[], 'TRB':[], 'TRG':[], 'IGH':[], 'IGL':[], 'IGK':[]}
    for idx, haplotype_data in enumerate(list_data):
        for data in list_data[idx]:
            allele_name = data[3]
            if 'OR' in allele_name:
                if 'TR' in allele_name:
                    dict_chain[idx]['TR/OR'].append(data)
                else:
                    dict_chain[idx]['IG/OR'].append(data)
            elif 'TRA' in allele_name or 'TRD' in allele_name:
                dict_chain[idx]['TRA/D'].append(data)
            elif 'TRB' in allele_name:
                dict_chain[idx]['TRB'].append(data)
            elif 'TRG' in allele_name:
                dict_chain[idx]['TRG'].append(data)
            elif 'IGH' in allele_name:
                dict_chain[idx]['IGH'].append(data)
            elif 'IGL' in allele_name:
                dict_chain[idx]['IGL'].append(data)
            elif 'IGK' in allele_name:
                dict_chain[idx]['IGK'].append(data)
            else:
                print('Error Name:', allele_name)

    # output according to gene chain
    if fn_out == None:  # print out the result
        for idx in range(2):
            print("==================================haplotype " + str(idx+1) + "========================================")
            for gene_chain, list_data in dict_chain[idx].items():
                print(gene_chain)
                old_pos = -1
                if len(list_data) > 0:
                    old_ref = sorted(list_data)[0][0]
                for item in sorted(list_data):
                    ref_name, start, length, allele_name = item
                    if old_ref != ref_name:
                        print()
                
                    if old_pos == start and old_ref == ref_name:
                        print(ref_name + ' ' + str(start) + ' ' + str(start + length) + ' ' + allele_name + '\t*')
                    else:
                        print(ref_name + ' ' + str(start) + ' ' + str(start + length) + ' ' + allele_name)
                    old_pos = start
                    old_ref = ref_name
    else:   # write to two haplotypes bed files
        list_out_name = [fn_out+".1.bed", fn_out+".2.bed"]
        for idx in range(2):
            output_string = ""
            for gene_chain, list_data in dict_chain[idx].items():
                output_string += (gene_chain + '\n')
                if len(list_data) > 0:
                    old_ref = sorted(list_data)[0][0]
                    old_name = []
                    old_pos = -1
                    old_len = -1
                    for item in sorted(list_data):
                        ref_name, start, length, allele_name = item
                        if old_ref == ref_name:
                            if old_pos + old_len >= start:
                                old_name.append(allele_name)
                                old_len = max(old_len, length)
                            else:
                                if old_pos != -1:
                                    output_string += (ref_name + ' ' + str(old_pos) + ' ' + str(old_pos + old_len) + ' ' + ';'.join(sorted(old_name)) + '\n')
                                
                                old_name = [allele_name]
                                old_pos = start
                                old_len = length
                        else:
                            if old_pos != -1:
                                output_string += (old_ref + ' ' + str(old_pos) + ' ' + str(old_pos + old_len) + ' ' + ';'.join(sorted(old_name)) + '\n')
                                output_string += ('\n')
                            old_name = [allele_name]
                            old_pos = start
                            old_len = length
                        old_ref = ref_name
                    if old_pos != -1:
                        output_string += (old_ref + ' ' + str(old_pos) + ' ' + str(old_pos + old_len) + ' ' + ';'.join(sorted(old_name)) + '\n')
            f = open(list_out_name[idx], 'w')
            f.write(output_string)
            f.close()



if __name__ == "__main__":
    main()

