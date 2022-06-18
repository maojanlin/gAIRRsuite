import pysam
import argparse

TRA_chain = {'TRA', 'TRD'}
TRB_chain = {'TRB'}
TRG_chain = {'TRG'}
IGH_chain = {'IGH'}
IGL_chain = {'IGL'}
IGK_chain = {'IGK'}





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--list_report', nargs='+', default=[], help='list of the report files of all the TR or IG')
    parser.add_argument('-o', '--output', help='the output report file')
    args = parser.parse_args()
    
    list_report = args.list_report
    fn_out      = args.output

    # read all the data files
    list_data = []
    for report_name in list_report:
        f = open(report_name, 'r')
        for line in f:
            fields = line.strip().split(',')
            if len(fields) < 2:
                continue
            elif len(fields) > 5:
                allele_name = fields[0] + '(' + fields[4][3:] + ',h' + fields[5][4:] + ')'
            elif len(fields) > 4:
                allele_name = fields[0] + '(' + fields[4][3:] + ')'
            else:
                allele_name = fields[0]
            list_data.append((fields[1],int(fields[2]),int(fields[3]),allele_name))

    # sort the data files
    dict_chain = {'TR/OR':[], 'IG/OR':[], 'TRA/D':[], 'TRB':[], 'TRG':[], 'IGH':[], 'IGL':[], 'IGK':[]}
    for data in list_data:
        allele_name = data[3]
        if 'OR' in allele_name:
            if 'TR' in allele_name:
                dict_chain['TR/OR'].append(data)
            else:
                dict_chain['IG/OR'].append(data)
        elif 'TRA' in allele_name or 'TRD' in allele_name:
            dict_chain['TRA/D'].append(data)
        elif 'TRB' in allele_name:
            dict_chain['TRB'].append(data)
        elif 'TRG' in allele_name:
            dict_chain['TRG'].append(data)
        elif 'IGH' in allele_name:
            dict_chain['IGH'].append(data)
        elif 'IGL' in allele_name:
            dict_chain['IGL'].append(data)
        elif 'IGK' in allele_name:
            dict_chain['IGK'].append(data)
        else:
            print('Error Name:', allele_name)

    # output according to gene chain
    for gene_chain, list_data in dict_chain.items():
        print(gene_chain)
        old_pos = -1
        for item in sorted(list_data):
            ref_name, start, length, allele_name = item
            if old_pos == start:
                print(ref_name + ' ' + str(start) + ' ' + str(start + length) + ' ' + allele_name + '\t*')
            else:
                print(ref_name + ' ' + str(start) + ' ' + str(start + length) + ' ' + allele_name)
            old_pos = start




