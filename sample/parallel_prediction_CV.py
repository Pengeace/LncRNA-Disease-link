import multiprocessing
from functools import reduce
from math import sqrt

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

fold = 5
calcu_sim = 0
normalize_score = 0
lncRNA_disease_path = '../data/lncRNA-disease.csv'
reserved_lncRNA_disease_pairs_path = '../data/data-partition/moved_pairs.txt'
item_index_path = '../data/item_index.txt'
data_partition_dir = '../data/data-partition/'


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


def load_train_test_split(path):
    train_test = [[], []]
    f = open(path, 'r')
    flag = 0
    for line in f.readlines():
        line = line.strip()
        if 'Train' in line:
            flag = 0
        elif 'Test' in line:
            flag = 1
        else:
            lncRNA, disease, label = [int(x) for x in line.split()]
            train_test[flag].append([lncRNA, disease, label])
    f.close()
    return train_test[0], train_test[1]


def load_reserved_pairs(path):
    reserved_pairs = []
    with open(path, 'r') as f:
        for line in f:
            reserved_pairs.append([int(x) for x in line.split()])
    return reserved_pairs


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


def process_each_fold(cur_fold, auc_fold, auc_disease_fold, lock):
    print("# Fold %d" % cur_fold)
    train, test = load_train_test_split(data_partition_dir + 'partition_fold{}.txt'.format(cur_fold))
    embeddings, node_num, dim = load_embedding(data_partition_dir + "embeddings_fold{}.txt".format(cur_fold))

    disease_related_lncRNAs = {}
    for [lncRNA, disease, label] in train + reserved_pairs:
        # for [lncRNA, disease, label] in train:
        d = index2name[disease]
        if d not in disease_related_lncRNAs:
            disease_related_lncRNAs[d] = []
        else:
            disease_related_lncRNAs[d].append([index2name[lncRNA], label])

    if calcu_sim:
        ll_sim = pd.DataFrame(np.zeros([len(lncRNAs), len(lncRNAs)]), index=lncRNAs, columns=lncRNAs)
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
        ll_sim.to_csv(data_partition_dir + 'lncRNA_sim{}.txt'.format(cur_fold))
    else:
        ll_sim = pd.read_csv(data_partition_dir + 'lncRNA_sim{}.txt'.format(cur_fold)).set_index('Unnamed: 0')

    ld = pd.read_csv(lncRNA_disease_path)
    ld = ld.set_index("Name")

    probas = []
    labels = []
    probas_disease = {}
    labels_disease = {}
    print("Predicting...")
    for pair in test:
        lncRNA, disease, label = pair
        lncRNA = index2name[lncRNA]
        disease = index2name[disease]
        sims = list(ll_sim.ix[lncRNA, [x[0] for x in disease_related_lncRNAs[disease]]])
        sum_sims = sum(sims)
        labs = list([x[1] for x in disease_related_lncRNAs[disease]])
        pred = (vector_multiply(sims, labs) / sum_sims) if sum_sims > 0 else 0.0

        probas.append(pred)
        labels.append(label)

        if disease not in probas_disease:
            probas_disease[disease] = []
            labels_disease[disease] = []
        probas_disease[disease].append(pred)
        labels_disease[disease].append(label)

    lock.acquire()
    auc_fold.append(roc_auc_score(labels, probas))
    print(auc_fold)
    for d in probas_disease:
        if max(labels_disease[d]) > 0:
            if d not in auc_disease_fold:
                auc_disease_fold[d] = []
            auc_disease_fold[d] = auc_disease_fold[d] + [roc_auc_score(labels_disease[d], probas_disease[d])]
    lock.release()


index2name, name2index, item_list = read_index_file(item_index_path)
diseases = [x[1] for x in item_list if x[2] == 'disease']
lncRNAs = [x[1] for x in item_list if x[2] == 'lncRNA']
reserved_pairs = load_reserved_pairs(reserved_lncRNA_disease_pairs_path)

if __name__ == '__main__':

    auc_fold = multiprocessing.Manager().list()
    auc_disease_fold = multiprocessing.Manager().dict()
    lock = multiprocessing.Lock()
    plist = [multiprocessing.Process(target=process_each_fold, args=[f, auc_fold, auc_disease_fold, lock]) for f in
             range(fold)]

    for p in plist:
        p.start()
    for p in plist:
        p.join()

    print("# AUC for all diseases:")
    print(auc_disease_fold)
    print("AUC fold CV:", np.mean(auc_fold))
    print("# AUC for all diseases in CV:")
    for d in auc_disease_fold.keys():
        print(d, np.mean(auc_disease_fold[d]))

# Previous:
# [0.9211684979426518, 0.9112173859429493, 0.9171714361079301, 0.9173682640220286, 0.9074345717172032]
# ('AUC fold CV:', 0.9148720311465525)
