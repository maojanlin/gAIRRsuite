
if __name__ == "__main__":
    # annotation ground truth
    set_annotation_H1 = set()
    f = open("flanking_region_analysis/annotation_NA12878_tcrv_contigs_H1.txt", "r")
    for line in f:
        set_annotation_H1.add(line.split(',')[1])
    print("Number of H1 alleles:", len(set_annotation_H1))

    set_annotation_H2 = set()
    f = open("flanking_region_analysis/annotation_NA12878_tcrv_contigs_H2.txt", "r")
    for line in f:
        set_annotation_H2.add(line.split(',')[1])
    print("Number of H2 alleles:", len(set_annotation_H2))

    # called flanking region
    set_remained_H1 = set()
    list_corrected_H1 = []
    f = open("target_call/NA12878_TCRV_flanking/flanking_result/bwa_flanking_H1.sam", "r")
    for line in f:
        if line[0] != '@':
            if line.split()[11] == "NM:i:0":
                list_corrected_H1.append(line.split()[0].split('/')[0])
            else:
                set_remained_H1.add(line.split()[0])
    set_corrected_H1 = set(list_corrected_H1)
    print("Number of correct H1 called:", len(set_corrected_H1))

    set_remained_H2 = set()
    list_corrected_H2 = []
    f = open("target_call/NA12878_TCRV_flanking/flanking_result/bwa_flanking_H2.sam", "r")
    for line in f:
        if line[0] != '@':
            if line.split()[11] == "NM:i:0":
                list_corrected_H2.append(line.split()[0].split('/')[0])
            else:
                set_remained_H2.add(line.split()[0])
    set_corrected_H2 = set(list_corrected_H2)
    print("Number of correct H2 called:", len(set_corrected_H2))

    print("FN in H1:", set_annotation_H1 - set_corrected_H1 )
    print("FP in H1:", set_corrected_H1 - set_annotation_H1 )
    print("FN in H2:", set_annotation_H2 - set_corrected_H2 )
    print("FP in H2:", set_corrected_H2 - set_annotation_H2 )
    print("There are", len(set_remained_H1.intersection(set_remained_H2)) ,"useless seq:")
    for element in sorted(set_remained_H1.intersection(set_remained_H2)):
        print('\t',element)


    


