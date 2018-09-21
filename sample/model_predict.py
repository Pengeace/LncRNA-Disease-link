import pandas as pd
import numpy as np
from math import sqrt
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold


fold = 5
lncRNA_disease_path = '../data/lncRNA-disease.csv'
item_index_path = '../data/item_index.txt'
data_partition_dir = '../data/data-partition/'


def vector_multiply(A, B):
    my_sum = lambda x, y: x + y
    return reduce(my_sum, [a*b for a,b in zip(A, B)])

def cosine_similarity(A,B):
    t1 = vector_multiply(A, B)
    t2 = sqrt(vector_multiply(A, A))
    t3 = sqrt(vector_multiply(B, B))
    if t2>0 and t3>0:
        return t1/(t2*t3)
    else:
        return 0


def load_train_test_split(path):
    train_test = [[],[]]
    f = open(path,'r')
    flag = 0
    for line in f.readlines():
        line = line.strip()
        if 'Train' in line:
            flag = 0
        elif 'Test' in line:
            flag = 1
        else:
            lncRNA, disease, label = [ int(x) for x in line.split()]
            train_test[flag].append([lncRNA, disease, label])
    return train_test[0], train_test[1]

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
diseases = [x[1] for x in item_list if x[2]=='disease']
lncRNAs = [x[1] for x in item_list if x[2]=='lncRNA']


# k-fold cross validation

auc_fold = []
auc_disease_fold = {}
for cur_fold in range(fold):
    print("# Fold %d" % cur_fold)
    train, test = load_train_test_split(data_partition_dir+'partition_fold{}.txt'.format(cur_fold))
    embeddings, node_num, dim = load_embedding(data_partition_dir + "embeddings_fold{}.txt".format(cur_fold))

    ll_sim = pd.DataFrame(np.zeros([len(lncRNAs), len(lncRNAs)]), index=lncRNAs, columns=lncRNAs)

    ld = pd.read_csv(lncRNA_disease_path)
    ld = ld.set_index("Name")

    for [lncRNA, disease, label] in test:
        if label==1:
            ld[index2name[lncRNA], index2name[disease]] = 0

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

    sum_lncRNA_sim = {}
    for lncRNA in lncRNAs:
        sum_lncRNA_sim[lncRNA] = sum(list(ll_sim[lncRNA]))


    probas = []
    labels = []
    probas_disease = {}
    labels_disease = {}

    for pair in test:
        lncRNA, disease, label = pair
        lncRNA = index2name[lncRNA]
        disease = index2name[disease]
        pred = (vector_multiply(list(ll_sim[lncRNA]), list(ld.ix[:,disease])) / sum_lncRNA_sim[lncRNA]) if sum_lncRNA_sim[lncRNA]>0 else 0.0

        probas.append(pred)
        labels.append(label)

        if disease not in probas_disease:
            probas_disease[disease] = []
            labels_disease[disease] = []
        probas_disease[disease].append(pred)
        labels_disease[disease].append(label)

    auc_fold.append(roc_auc_score(labels, probas))

    print(auc_fold)
    for d in probas_disease:
        if max(labels_disease[d]) > 0:
            if d not in auc_disease_fold:
                auc_disease_fold[d] = []
            auc_disease_fold[d].append(roc_auc_score(labels_disease[d], probas_disease[d]))

print("# AUC for all diseases:")
print(auc_disease_fold)
print("AUC fold CV:",np.mean(auc_fold))
print("# AUC for all diseases in CV:")
for d in auc_disease_fold:
    print(d,np.mean(auc_disease_fold[d]))

# AUC in 5 fold
# [0.9154917976976293, 0.9186329353080193, 0.912639027764002, 0.9112176805252827, 0.9284490649617969]
# Mean : 0.9172861012513461
