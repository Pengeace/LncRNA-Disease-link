import re
import pandas as pd

lncRNA_disease_path = '../data/lncRNA_disease.csv'
microRNA_disease_path = '../data/microRNA_disease.csv'
microRNA_lncRNA_path = '../data/microRNA_lncRNA.csv'
item_index_path = '../data/item_index.txt'


def read_index_file(path):
    index2name = {}
    name2index = {}
    item_list = []
    with open(path) as f:
        for line in f:
            item = line.strip().split('\t')
            item_list.append(item)
            index2name[int(item[0])]= item[1]
            name2index[item[1]] = int(item[0])
    return index2name, name2index, item_list

def find_duplicated_items(target_list, name_map):
    duplicated_pairs = []
    for i in range(len(target_list)):
        x = target_list[i]
        for j in range(i + 1, len(target_list)):
            y = target_list[j]
            if name_map[x] == name_map[y]:
                duplicated_pairs.append([x, y])
    return duplicated_pairs

index2name, name2index, item_list = read_index_file(item_index_path)

lncRNAs = []
microRNAs = []
diseases = []
name_map = {}

for it in item_list:
    if it[2]=='lncRNA':
        lncRNAs.append(it[1])
    elif it[2]=='microRNA':
        microRNAs.append(it[1])
    elif it[2]=='disease':
        diseases.append(it[1])
    name_map[it[1]] = ''.join([x for x in it[1].lower() if x.isalnum()])

duplicated_lncRNAs = find_duplicated_items(lncRNAs, name_map)
duplicated_microRNAs = find_duplicated_items(microRNAs, name_map)
duplicated_diseases = find_duplicated_items(diseases, name_map)


for x in lncRNAs:
    if '(' in x and ')' in x:
        y = re.findall('[(](.*)[)]',x)[0]
        if y in lncRNAs:
            duplicated_lncRNAs.append([x, y])

print('# Duplicated lncRNAs pairs (%d):' % (len(duplicated_lncRNAs)))
for x in duplicated_lncRNAs:
    print(x)

print('# Duplicated microRNAs (%d):' % (len(duplicated_microRNAs)) )
for x in duplicated_microRNAs:
    print(x)

print('# Duplicated diseases (%d):' % (len(duplicated_diseases)))
for x in duplicated_diseases:
    print(x)