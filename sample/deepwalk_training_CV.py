import os
import random
from copy import deepcopy

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

fold = 5
random_seed = 0
data_partition_dir = '../data/data-partition/'
lncRNA_disease_path = '../data/lncRNA-disease.csv'
microRNA_disease_path = '../data/microRNA-disease.csv'
microRNA_lncRNA_path = '../data/microRNA-lncRNA.csv'
item_index_path = '../data/item_index.txt'


def read_index_file(path):
    index2name = {}
    name2index = {}
    item_list = []
    with open(path) as f:
        for line in f:
            item = line.strip().split('\t')
            item_list.append(item)
            index2name[int(item[0])] = item[1]
            if item[1] not in name2index:
                name2index[item[1]] = int(item[0])
            else:
                print(item[1])
    return index2name, name2index, item_list


def update_adjacent_list(data_frame, name2index, adj_list):
    columns = list(data_frame.columns)[1:]
    rows = list(data_frame['Name'])
    data_frame = data_frame.set_index('Name')

    for r in rows:
        item = name2index[r]
        adjs = dict(data_frame.ix[r])
        if item not in adj_list:
            adj_list[item] = []
        for key in adjs:
            if adjs[key]:
                adj_list[item].append(name2index[key])

    for c in columns:
        item = name2index[c]
        adjs = dict(data_frame[c])
        if item not in adj_list:
            adj_list[item] = []
        for key in adjs:
            if adjs[key]:
                adj_list[item].append(name2index[key])

    return adj_list


index2name, name2index, item_list = read_index_file(item_index_path)
diseases = [x[1] for x in item_list if x[2] == 'disease']
lncRNAs = [x[1] for x in item_list if x[2] == 'lncRNA']
microRNA_indexes = [x[0] for x in item_list if x[2] == 'microRNA']
ld = pd.read_csv(lncRNA_disease_path)
md = pd.read_csv(microRNA_disease_path)
ml = pd.read_csv(microRNA_lncRNA_path)

base_adj_list = {}
for data in [md, ml]:
    base_adj_list = update_adjacent_list(data, name2index, base_adj_list)

ld_pos_adj_list = update_adjacent_list(ld, name2index, {})

total_pairs = []
positive_pairs = []
ld = ld.set_index('Name')
for l in ld.index:
    for d in ld.columns:
        label = int(ld.ix[l, d])
        total_pairs.append([name2index[l], name2index[d], label])
        if label == 1:
            positive_pairs.append([name2index[l], name2index[d], label])


# Create the minimal linked network
#############################

def find(x, parents):
    if parents[x] < 0:
        return x
    else:
        p = find(parents[x], parents)
        parents[x] = p
        return p


parents = [-1] * len(name2index.values())
for key in base_adj_list:
    for val in base_adj_list[key]:
        x = find(key, parents)
        y = find(val, parents)
        if x != y:
            if x in microRNA_indexes:
                parents[x] = y
            else:
                parents[y] = x

moved_pairs = []
clusters = [x for x in range(len(parents)) if parents[x] < 0]
random.seed(random_seed)
print("Clusters:", clusters)
while len(clusters) > 1:
    random_pair = positive_pairs[random.randint(0, len(positive_pairs) - 1)]
    x = find(random_pair[0], parents)
    y = find(random_pair[1], parents)
    if x != y:
        # print("Dealing with {}, cluster size : {}".format(random_pair, len(clusters)))

        moved_pairs.append(random_pair)
        total_pairs.remove(random_pair)

        base_adj_list[random_pair[0]].append(random_pair[1])
        base_adj_list[random_pair[1]].append(random_pair[0])

        ld_pos_adj_list[random_pair[0]].remove(random_pair[1])
        ld_pos_adj_list[random_pair[1]].remove(random_pair[0])

        parents[x] = y

    clusters = [x for x in range(len(parents)) if parents[x] < 0]

with open(data_partition_dir + 'moved_pairs.txt', 'w') as f:
    for item in moved_pairs:
        f.write('\t'.join([str(x) for x in item]) + '\n')

##############################

# k-fold cross validation
skf = StratifiedKFold(n_splits=fold, shuffle=True, random_state=random_seed)
cur_fold = 0
total_pairs = np.array(total_pairs)
for train, test in skf.split(total_pairs, [x[2] for x in total_pairs]):
    print('# fold %d' % cur_fold)
    pairs_train = total_pairs[train]
    pairs_test = total_pairs[test]

    print("# Writing data-partition pattern...")
    with open(data_partition_dir + 'partition_fold{}.txt'.format(cur_fold), 'w') as f:
        f.write("# Train set\t%d\n" % (len(pairs_train)))
        for item in pairs_train:
            f.write('\t'.join([str(x) for x in item]) + '\n')
        f.write('# Test set\t%d\n' % (len(pairs_test)))
        for item in pairs_test:
            f.write('\t'.join([str(x) for x in item]) + '\n')

    adj_list = deepcopy(base_adj_list)
    for [lncRNA, disease, label] in pairs_train:
        if label:
            adj_list[lncRNA].append(disease)
            adj_list[disease].append(lncRNA)

    with open(data_partition_dir + '/MLD_network_fold{}.txt'.format(cur_fold), 'w') as f:
        for key in adj_list:
            f.write('\t'.join([str(key)] + [str(x) for x in adj_list[key]]) + '\n')

    print("Deepwalk training...")
    os.system("deepwalk --input " + data_partition_dir + '/MLD_network_fold{}.txt '.format(cur_fold)
              + "--number-walks 80 --representation-size 128 "
              + "--walk-length 40 --window-size 10 --workers 12 --output " + data_partition_dir + '/embeddings_fold{}.txt'.format(
        cur_fold))
    cur_fold = cur_fold + 1

os.system('python ./parallel_prediction_CV.py')
