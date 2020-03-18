gene = 'TRBV'

list_blastn = []
# with open('trbv_top50.txt', 'r') as f:
with open('trbv_top80.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if line.count(gene) == 0:
            print ('warning:', line)
        else:
            list_blastn.append(line)

list_annotated = []
with open('trbv_all.txt', 'r') as f:
    for line in f:
        line = line.rstrip()
        if line.count(gene) == 0:
            print ('warning:', line)
        else:
            list_annotated.append(line)

set_blastn = set(list_blastn)
set_annotated = set(list_annotated)

print ('Size of blastn alleles:', len(set_blastn))
print ('Size of annotated alleles:', len(set_annotated))
print ('Size of intersection:', len(set_blastn.intersection(set_annotated)))
print ('Size of union:', len(set_blastn.union(set_annotated)))

print ('Found by blastn, not annotated:', set_blastn - set_blastn.intersection(set_annotated))
