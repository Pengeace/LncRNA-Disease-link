import os
import pandas



ld = pandas.read_csv('./data-union/total_lncRNADisease2.csv')   # 1831 x 281




lncRNAs = list(ld['Name'])  # 1831
diseases = list(ld.columns)[1:] # 281

print(len(lncRNAs),len(diseases))

with open('./data-union/item_index.txt', 'w') as f:
    index = 1
    for lncRNA in lncRNAs:
        f.write('%d\t%s\t%s\n' % (index, lncRNA, 'lncRNA'))
        index = index + 1
    for disease in diseases:
        f.write('%d\t%s\t%s\n' % (index, disease, 'disease'))
        index = index + 1