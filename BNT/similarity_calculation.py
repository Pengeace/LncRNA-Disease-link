import pandas
from functools import reduce
from math import sqrt
from sklearn.metrics import roc_curve, auc
import numpy as np
import os

def load_embedding(path):
    head = False
    num = None
    dim = None
    embeddings = {}
    with open(path,'r') as f:
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
    return reduce(my_sum ,[a*b for a,b in zip(A, B)])

def cosine_similarity(A,B):
    t1 = vector_multiply(A, B)
    t2 = sqrt(vector_multiply(A, A))
    t3 = sqrt(vector_multiply(B, B))
    if t2>0 and t3>0:
        return t1/(t2*t3)
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

index2name, name2index, item_list = read_index_file('./data-union/item_index.txt')
diseases = [x[1] for x in item_list if x[2]=='disease']
lncRNAs = [x[1] for x in item_list if x[2]=='lncRNA']

embeddings, num, dim = load_embedding('./data-union/embedding.txt')
ld = pandas.read_csv('./data-union/total_lncRNADisease2.csv')
ld = ld.set_index("Name")   # 2422 x 281


if not os.path.exists('disease_disease_similarity.csv'):

    dd_sim = pandas.DataFrame(np.zeros([len(diseases), len(diseases)]), index=diseases, columns=diseases)
    ll_sim = pandas.DataFrame(np.zeros([len(lncRNAs), len(lncRNAs)]), index=lncRNAs, columns=lncRNAs)

    # calculate disease-disease similarity
    print("Calculating disease-disease similarities...")
    for d1 in diseases:
        for d2 in diseases:
            if d1==d2:
                dd_sim.ix[d1,d2] = 0
            else:
                sim = cosine_similarity(embeddings[name2index[d1]], embeddings[name2index[d2]])
                dd_sim.ix[d1, d2] = sim
                dd_sim.ix[d2, d1] = sim


    # calculate lncRNA-lncRNA similarity
    print("Calculating lncRNA-lncRNA similarities...")
    for l1 in lncRNAs:
        for l2 in lncRNAs:
            if l1==l2:
                ll_sim.ix[l1,l2] = 0
            else:
                sim = cosine_similarity(embeddings[name2index[l1]], embeddings[name2index[l2]])
                ll_sim.ix[l1, l2] = sim
                ll_sim.ix[l2, l1] = sim

    dd_sim.to_csv('disease_disease_similarity.csv')
    ll_sim.to_csv('lncRNA_lncRNA_similarity.csv')

else:
    dd_sim = pandas.read_csv('disease_disease_similarity.csv')
    ll_sim = pandas.read_csv('lncRNA_lncRNA_similarity.csv')



lncRNA_disease_link_DBSI = pandas.DataFrame(np.zeros([len(lncRNAs),len(diseases)]), index=lncRNAs, columns=diseases)
lncRNA_disease_link_TBSI = pandas.DataFrame(np.zeros([len(lncRNAs),len(diseases)]), index=lncRNAs, columns=diseases)

# DBSI process
print("# DBSI process...")
y_labels = []
y_probas = []

for disease in diseases:
    for lncRNA in lncRNAs:
        lncRNA_disease_link_DBSI.ix[lncRNA, disease] = vector_multiply(list(dd_sim[disease]), list(ld.ix[lncRNA])) / vector_multiply(list(dd_sim[disease]), [1]*len(diseases))
    before_norm = np.array(lncRNA_disease_link_DBSI.ix[:,disease])
    lncRNA_disease_link_DBSI.ix[:, disease] = (before_norm - before_norm.min()) / before_norm.max() - before_norm.min()


for lncRNA in lncRNAs:
    for disease in diseases:
        y_labels.append(ld.ix[lncRNA,disease])
        y_probas.append(lncRNA_disease_link_DBSI.ix[lncRNA,disease])
fpr, tpr, thresholds = roc_curve(y_labels, y_probas)
AUC = auc(fpr, tpr)
print("DBSI AUC value: %f" % AUC)


# TBSI process
print("# TBSI process...")
y_labels = []
y_probas = []

for lncRNA in lncRNAs:
    for disease in diseases:
        lncRNA_disease_link_TBSI.ix[lncRNA, disease] = vector_multiply(list(ll_sim[lncRNA]), list(ld.ix[:,disease])) / vector_multiply(list(ll_sim[lncRNA]), [1]*len(lncRNAs))
    before_norm = np.array(lncRNA_disease_link_TBSI.ix[lncRNA,:])
    lncRNA_disease_link_TBSI.ix[lncRNA,:] = (before_norm - before_norm.min()) / before_norm.max() - before_norm.min()


for lncRNA in lncRNAs:
    for disease in diseases:
        y_labels.append(ld.ix[lncRNA,disease])
        y_probas.append(lncRNA_disease_link_TBSI.ix[lncRNA,disease])
fpr, tpr, thresholds = roc_curve(y_labels, y_probas)
AUC = auc(fpr, tpr)
print("TBSI AUC value: %f" % AUC)


# Calculating disease-disease similarities...
# Calculating lncRNA-lncRNA similarities...
# # DBSI process...
# DBSI AUC value: 0.501467
# # TBSI process...
# TBSI AUC value: 0.911197