import pandas

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
            index2name[int(item[0])]= item[1]
            name2index[item[1]] = int(item[0])

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

def get_sub_network(start, net, flags, adj_list):
    net.append(start)
    flags[start] = 0
    for x in adj_list[start]:
        if flags[x] > 0:
            net = get_sub_network(x, net, flags, adj_list)
    return net

ld = pandas.read_csv(lncRNA_disease_path)
md = pandas.read_csv(microRNA_disease_path)
ml = pandas.read_csv(microRNA_lncRNA_path)


microRNAs = list(md['Name'])
lncRNAs = list(ld['Name'])
diseases = list(ld.columns[1:])

index2name, name2index, item_list = read_index_file(item_index_path)
adj_list = {}
for data in [ld, md, ml]:
    adj_list = update_adjacent_list(data, name2index, adj_list)


flags = [1] * len(item_list)    # whether a item has been included into a network
networks = []
for i in range(len(flags)):
    if flags[i]>0:
        # print(sum(flags),'...')
        networks.append(get_sub_network(i,[],flags,adj_list))

num_network = len(networks)
biggest_network = None
print(num_network)
for i in range(num_network):
    if len(networks[i])>2000:
        biggest_network = networks[i]
    # print('# Network %d, size %d:' % (i+1, len(networks[i])))
    for x in networks[i]:
        print(item_list[x])

# lncRNAs_new = []
# microRNAs_new = []
# diseases_new = []
#
#
# for x in biggest_network:
#     item = item_list[x]
#     if item[2]=='microRNA':
#         microRNAs_new.append(item[1])
#     if item[2]=='lncRNA':
#         lncRNAs_new.append(item[1])
#     if item[2]=='disease':
#         diseases_new.append(item[1])
#
# print(len(microRNAs_new), len(lncRNAs_new), len(diseases_new))
# # (675, 1723, 236)
# lncRNAs_removed = [x for x in lncRNAs if x not in lncRNAs_new]
# microRNAs_removed = [x for x in microRNAs if x not in microRNAs_new]
# diseases_removed = [x for x in diseases if x not in diseases_new]
#
#
# ld = ld.set_index('Name')
# md = md.set_index('Name')
# ml = ml.set_index('Name')
#
# ld = ld.drop(lncRNAs_removed, axis=0)
# ld = ld.drop(diseases_removed, axis=1)
# md = md.drop(microRNAs_removed, axis=0)
# md = md.drop(diseases_removed, axis=1)
# ml = ml.drop(microRNAs_removed, axis=0)
# ml = ml.drop(lncRNAs_removed, axis=1)
#
# lncRNA_disease_path = '../data/lncRNA-disease.csv'
# microRNA_disease_path = '../data/microRNA-disease.csv'
# microRNA_lncRNA_path = '../data/microRNA-lncRNA.csv'
#
# ld.to_csv(lncRNA_disease_path)
# md.to_csv(microRNA_disease_path)
# ml.to_csv(microRNA_lncRNA_path)