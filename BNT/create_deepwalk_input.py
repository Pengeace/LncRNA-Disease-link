import pandas


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

def get_adjacent_list(data_frame, name2index):
    columns = list(data_frame.columns)[1:]
    rows = list(data_frame['Name'])
    data_frame = data_frame.set_index('Name')
    adj_list = []
    for r in rows:
        item = [name2index[r]]
        adjs = dict(data_frame.ix[r])
        for key in adjs:
            if adjs[key]:
                item.append(name2index[key])
        if len(item)>1:
            adj_list.append(item)

    for c in columns:
        item = [name2index[c]]
        adjs = dict(data_frame[c])
        for key in adjs:
            if adjs[key]:
                item.append(name2index[key])
        if len(item)>1:
            adj_list.append(item)
    return adj_list


ld = pandas.read_csv('./data-union/total_lncRNADisease2.csv')   # 2433 x 282


index2name, name2index, item_list = read_index_file('./data-union/item_index.txt')

adj_list = []
for data in [ld]:
    adj_list = adj_list + get_adjacent_list(data, name2index)

with open('./data-union/lncRNA-disease_network.txt','w') as f:
    for item in adj_list:
        f.write('\t'.join([str(x) for x in item])+'\n')