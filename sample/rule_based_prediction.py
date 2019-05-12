import os
from functools import reduce
from math import sqrt

import numpy as np
import pandas
from sklearn.metrics import roc_auc_score

data_base_dir = '../data/'
item_index_path = data_base_dir + 'item_index.txt'
lncRNA_disease_path = data_base_dir + 'lncRNA-disease.csv'
lncRNA_sim_path = data_base_dir + 'lncRNA_lncRNA_similarity.csv'
lncRNA_disease_prediction_path = data_base_dir + 'lncRNA-disease-prediction.csv'
nomalize_score = False


def load_embedding(path):
    head = False
    num = None
    dim = None
    embeddings = {}
    with open(path, 'r') as f:
        for line in f:
            if head:
                num, dim = [int(x) for x in line.strip().split()]
                head = False
            else:
                items = line.strip().split()
                embeddings[int(items[0])] = [float(x) for x in items[1:]]
    return embeddings, num, dim


def vector_multiply(A, B):
    my_sum = lambda x, y: x + y
    return reduce(my_sum, [a * b for a, b in zip(A, B)])


def cosine_similarity(A, B):
    t1 = vector_multiply(A, B)
    t2 = sqrt(vector_multiply(A, A))
    t3 = sqrt(vector_multiply(B, B))
    if t2 > 0 and t3 > 0:
        return t1 / (t2 * t3)
    else:
        return 0


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


index2name, name2index, item_list = read_index_file(item_index_path)
diseases = [x[1] for x in item_list if x[2] == 'disease']
lncRNAs = [x[1] for x in item_list if x[2] == 'lncRNA']

ld = pandas.read_csv(lncRNA_disease_path)
ld = ld.set_index("Name")

if not os.path.exists(lncRNA_sim_path):

    embeddings, num, dim = load_embedding('../data/embeddings.txt')

    ll_sim = pandas.DataFrame(np.zeros([len(lncRNAs), len(lncRNAs)]), index=lncRNAs, columns=lncRNAs)

    # calculate lncRNA-lncRNA similarity
    print("Calculating lncRNA-lncRNA similarities...")
    for l1 in lncRNAs:
        for l2 in lncRNAs:
            if l1 == l2:
                ll_sim.ix[l1, l2] = 0
            else:
                sim = cosine_similarity(embeddings[name2index[l1]], embeddings[name2index[l2]])
                ll_sim.ix[l1, l2] = sim
                ll_sim.ix[l2, l1] = sim

    ll_sim.to_csv(lncRNA_sim_path)

else:
    ll_sim = pandas.read_csv(lncRNA_sim_path).set_index('Unnamed: 0')


lncRNA_disease_link_TBSI = pandas.DataFrame(np.zeros([len(lncRNAs), len(diseases)]), index=lncRNAs, columns=diseases)

# TBSI process
print("# Predict lncRNA_disease link based on lncRNA similarities...")
y_labels = []
y_probas = []

for lncRNA in lncRNAs:
    sum_lncRNA_sim = sum(list(ll_sim[lncRNA]))
    for disease in diseases:
        lncRNA_disease_link_TBSI.ix[lncRNA, disease] = (vector_multiply(list(ll_sim[lncRNA]), list(
            ld.ix[:, disease])) / sum_lncRNA_sim) if sum_lncRNA_sim > 0 else 0.0
    if nomalize_score:
        before_norm = np.array(lncRNA_disease_link_TBSI.ix[lncRNA, :])
        max_v = before_norm.max()
        min_v = before_norm.min()
        if max_v - min_v > 0:
            lncRNA_disease_link_TBSI.ix[lncRNA, :] = (before_norm - min_v) / (max_v - min_v)

for lncRNA in lncRNAs:
    for disease in diseases:
        y_labels.append(ld.ix[lncRNA, disease])
        y_probas.append(lncRNA_disease_link_TBSI.ix[lncRNA, disease])

# store result
lncRNA_disease_link_TBSI.to_csv(lncRNA_disease_prediction_path)

AUC = roc_auc_score(y_labels, y_probas)
print("AUC value: %f" % AUC)

auc_dict = {}
for disease in diseases:
    labels = []
    probas = []
    for lncRNA in lncRNAs:
        labels.append(ld.ix[lncRNA, disease])
        probas.append(lncRNA_disease_link_TBSI.ix[lncRNA, disease])

    if max(labels) > 0:
        auc_dict[disease] = roc_auc_score(labels, probas)
auc_dict = sorted(auc_dict.items(), lambda x, y: cmp(x[1], y[1]), reverse=True)
with open('AUC for diseases.txt', 'w') as f:
    f.write('Disease\tAUC\n')
    for it in auc_dict:
        num = int(sum(list(ld[it[0]])))
        f.write(str(it[0]) + '\t' + str(num) + '\t' + str(it[1]) + '\n')
