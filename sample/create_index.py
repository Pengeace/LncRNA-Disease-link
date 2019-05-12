import pandas

lncRNA_disease_path = '../data/lncRNA-disease.csv'
microRNA_disease_path = '../data/microRNA-disease.csv'
microRNA_lncRNA_path = '../data/microRNA-lncRNA.csv'

ld = pandas.read_csv(lncRNA_disease_path)
md = pandas.read_csv(microRNA_disease_path)
ml = pandas.read_csv(microRNA_lncRNA_path)


microRNAs = list(md['Name'])
lncRNAs = list(ld['Name'])
diseases = list(ld.columns[1:])

print(len(microRNAs),len(lncRNAs),len(diseases))
# (675, 1723, 236)

with open('../data/item_index.txt', 'w') as f:
    index = 0
    for microRNA in microRNAs:
        f.write('%d\t%s\t%s\n' % (index, microRNA, 'microRNA'))
        index = index + 1
    for lncRNA in lncRNAs:
        f.write('%d\t%s\t%s\n' % (index, lncRNA, 'lncRNA'))
        index = index + 1
    for disease in diseases:
        f.write('%d\t%s\t%s\n' % (index, disease, 'disease'))
        index = index + 1