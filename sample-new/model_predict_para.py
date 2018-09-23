import multiprocessing
from math import sqrt

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

fold = 10
calcu_sim = 1
lncRNA_disease_path = '../data/lncRNA-disease.csv'
# reserved_lncRNA_disease_pairs_path = '../data/data-partition/moved_pairs.txt'
item_index_path = '../data/item_index.txt'
data_partition_dir = '../data/data-partition-new/'


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
    # for [lncRNA, disease, label] in train + reserved_pairs:
    for [lncRNA, disease, label] in train:
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
# reserved_pairs = load_reserved_pairs(reserved_lncRNA_disease_pairs_path)

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

# 5 fold
# [0.6810548015344323]
# [0.6810548015344323, 0.6847487847228895]
# [0.6810548015344323, 0.6847487847228895, 0.6511380500190863]
# [0.6810548015344323, 0.6847487847228895, 0.6511380500190863, 0.673129475092283]
# [0.6810548015344323, 0.6847487847228895, 0.6511380500190863, 0.673129475092283, 0.6579033447084337]
# # AUC for all diseases:
# {'Burkitts lymphoma': [0.10526315789473684, 0.061111111111111116, 0.48011363636363635], 'esophageal adenocancer': [0.9740634005763689, 0.8600583090379009], 'pancreas cancer': [0.4261012762453685, 0.099406528189911, 0.18106508875739646, 0.7644508670520231, 0.6251638269986894], 'nasopharyngeal cancer': [0.45189504373177847, 0.2749485596707819, 0.3091894060995184, 0.6622807017543859], 'Autism spectrum disorder': [0.5], 'panic disorder': [0.5], 'acute lymphoblastic leukemia': [0.09860248447204967, 0.46120689655172414, 0.2619718309859155], 'Intrauterine Growth Restriction': [0.5], 'B-cell neoplasms': [0.7624633431085045, 0.8773584905660378, 0.7313432835820896], "Huntington's disease": [0.21216617210682492, 0.6963402571711176, 0.4342105263157895, 0.6405797101449275], 'B-cell lymphoma': [0.8276923076923077, 0.0, 0.1899441340782123], 'dilated cardiomyopathy': [0.5480769230769231, 0.7823834196891192], 'rhabdomyosarcoma': [0.5], 'polymyositis': [0.5], 'Leishmania': [0.5], 'basal-like breast cancer': [0.5], 'promyelocytic leukemia': [0.5], 'pre-eclampsia': [0.5], 'atherosclerosis': [0.8593314763231198, 0.0, 0.9464285714285714, 0.8567251461988304], 'multiple myeloma': [0.5316091954022988, 0.6725239616613419, 0.7122186495176849, 0.9611111111111111, 0.8185378590078329], 'gastric cancer': [0.4264311270125224, 0.4889214046822743, 0.3828125, 0.35588740584115236, 0.3637829912023461], 'Obesity': [0.6288819875776397, 0.07210031347962381], 'lymphoma': [0.2571428571428571, 0.8260869565217391, 0.5325325325325325, 0.7038216560509554], 'Wiedemann-Beckwith syndrome': [0.5], 'oesophageal adenocancer': [0.5], 'adolescent idiopathic scoliosis': [0.10951008645533145, 0.16567164179104477], 'DiGeorge syndrome': [0.5], 'postmenopausal osteoporosis': [0.5], 'chronic lymphocytic leukemia': [0.024566473988439308, 0.5667613636363636, 0.6496913580246914], 'autoimmune disease': [0.584958217270195, 0.025423728813559365], 'myocardial infarction': [0.943952802359882, 0.855072463768116], 'Wilms tumor': [0.8861111111111111, 0.9715189873417721, 1.0, 0.8087774294670846, 0.9329608938547486], 'gestational trophoblastic diseases': [0.5], 'renal cell cancer': [0.2920845272206304, 0.4247366632687733, 0.31004531722054385, 0.21891327063740854, 0.5917142857142856], 'Crohn disease': [0.5], 'Duchenne muscular dystrophy': [0.7647058823529411, 0.06574923547400613], 'acute myocardial infarction': [0.2122905027932961, 0.1753424657534247, 0.17891373801916932, 0.04016620498614959], 'ovarian cancer': [0.7134986225895317, 0.7670623145400594, 0.37317784256559766, 0.5112612612612613, 0.7306748466257669], 'bone diseases': [0.5], 'schizophrenia': [0.497134670487106, 0.8303834808259588, 0.882768361581921, 0.9752321981424149, 0.9828080229226361], 'laryngeal cancer': [0.5], "Hodgkin's lymphoma": [0.5], 'frontotemporal lobar degeneration': [0.5], "Alzheimer's disease": [0.7528248587570622, 0.7125382262996942, 0.7357954545454546, 0.22544642857142858], 'myeloproliferative polycythaemia vera': [0.5], 'bipolar disorder': [0.9912023460410557, 0.9893899204244032], 'diabetic nephropathy': [0.5], 'basal cell cancer': [0.5], 'dyskeratosis congenita': [0.22379603399433423, 0.043296089385474856], 'psychiatric disease': [0.45868945868945865, 0.0028571428571428914], 'cervical cancer': [0.9223880597014926, 0.629296875, 0.3440233236151603, 0.6558067375886525, 0.45810397553516824], 'gallbladder cancer': [0.947075208913649, 0.5991620111731844, 0.4198717948717949], 'pheochromocytoma': [0.9310344827586207, 0.9527777777777777], 'HCV': [0.5], 'Prader-Willi syndrome': [0.016901408450704203, 0.06976744186046513, 0.5178571428571429, 0.6954022988505747], 'acute megakaryoblastic leukemia': [0.9537313432835821, 0.9971830985915493], 'Neurodevelopmental syndromes associated with the SOX2 locus': [0.5], 'velocardiofacial syndrome': [0.5], 'neural tube defects': [0.5], 'cardiomyopathy': [0.5], 'endometrial cancer': [0.6041033434650457, 0.017543859649122806, 0.676392572944297, 0.023668639053254448, 0.03230769230769226], 'type 2 diabetes': [0.0720461095100865, 0.48282265552460535, 0.9161490683229814], 'medulloblastoma': [0.5], 'adrenocortical cancers': [0.5], 'meningioma': [0.5], 'hepatocellular cancer': [0.34645286686103005, 0.4042308796783766, 0.3234989648033126, 0.3097604692846668, 0.4583554376657825], 'growth restriction': [0.5], 'Congenital hyperinsulinism': [0.5], 'testicular cancer': [0.03794037940379408, 0.3088235294117647, 0.6882022471910112, 0.002915451895043719], 'aortic aneurysm': [0.5], 'intracranial aneurism': [0.06784660766961653, 0.8022284122562674, 0.02949852507374634], 'cancer': [0.6221562809099901, 0.9195402298850576, 0.7043524416135881, 0.33461538461538465, 0.5182325182325183], 'malignant pleural mesothelioma': [0.9515669515669516, 0.3609467455621302, 0.009090909090909038, 0.6463768115942029, 0.4669540229885057], 'AIDS': [0.7970149253731343, 0.8138138138138138, 0.45592705167173253, 0.5953079178885631, 0.6767810026385224], 'trophoblastic tumor': [0.5], 'pseudohypoparathyroidism type Ib': [0.5], 'malignant melanoma': [0.5818713450292398, 0.7485714285714286, 0.7935103244837758], 'neurofibromatosis type 1': [0.5], 'acute promyelocytic leukemia': [0.5314285714285714, 0.6726190476190476], 'gastric adenocancer': [0.6119302949061662, 0.36477987421383645, 0.02077363896848136, 0.15298885511651467, 0.5093411996066863], 'hepatoma': [0.5], 'biparental complete hydatidiform moles': [0.5], 'gastrointestinal stromal tumor': [0.5], 'chronic myeloid leukemia': [0.5], 'breast cancer': [0.44486803519061585, 0.46805111821086265, 0.5013927576601671, 0.6045407458563535, 0.49722442505947667], 'amyotrophic lateral sclerosis': [0.5], 'Diabetes': [0.5647873392680515, 0.8025568181818181], 'psoriasis': [0.8961748633879781, 0.28125, 0.7654320987654322], 'fibrosarcoma': [0.5], 'Glaucoma': [0.5], 'depression': [0.5], 'colon cancer': [0.31378299120234604, 0.5340718361866399, 0.5314782982380748, 0.5632312671147625, 0.47717622080679406], 'Pituitary adenoma': [0.45538243626062325, 0.9269005847953217, 0.7605421686746987], 'oesophagus cancer': [0.5], 'chronic myeloproliferative disorders': [0.5], 'hematopoiesis': [0.5], 'small-cell lung cancer': [0.8944444444444445, 0.7574850299401198], 'Esophageal squamous cell cancer': [0.653030303030303, 0.5476866883116882, 0.36292613636363635, 0.5723415132924335, 0.40060422960725073], 'myopia': [0.5], 'membranous nephropathy': [0.5], 'cleft lip': [0.5], 'squamous cell cancer': [0.5739035087719297, 0.6319248826291081, 0.7122093023255813, 0.6779661016949152, 0.2078313253012048], 'coronary disease': [0.5], 'papillary thyroid cancer': [0.25785714285714284, 0.38168352601156064, 0.01869158878504673, 0.4757799671592775], 'prostate cancer': [0.47558728696453256, 0.8052173913043478, 0.4650866462793068, 0.5185185185185185, 0.3321220930232558], 'Inclusion body myositis': [0.5], 'osteosarcoma': [0.7469418960244647, 0.16515609264853975, 0.6753731343283581, 0.3850241545893719, 0.034340659340659385], 'glioblastoma': [0.6681681681681682, 0.3908045977011494, 0.6869230769230769, 0.5807142857142856, 0.5421511627906977], 'hyperhomocysteinemia': [0.5], 'sarcoma': [0.5], 'tongue cancer': [0.9384057971014492, 0.973134328358209], 'Beckwith-Wiedemann syndrome': [0.9266862170087976, 0.8870056497175142], 'periodontitis': [0.9658119658119658, 0.9827089337175793], 'lymphocytic leukemia': [0.5], 'aging': [0.0, 0.7121661721068249], 'melanoma': [0.268018018018018, 0.724570200573066, 0.6276041666666666, 0.8882521489971347, 0.877906976744186], 'hepatocelluar cancer': [0.7554296506137865, 0.5716292134831461, 0.6782352941176469, 0.6818181818181819, 0.756838905775076], 'bladder cancer': [0.5128205128205128, 0.5959302325581396, 0.2735943775100402, 0.568204365079365, 0.5064157399486741], 'drug abuse': [0.5], 'acute T-lymphocytic leukemia': [0.5], 'rheumatoid arthritis': [0.5], 'myelodysplastic syndrome': [0.049844236760124616, 0.44744744744744747], 'retinoblastoma': [0.5], 'multiple sclerosis': [0.5], 'coronary heart disease': [0.5], 'neuroblastoma': [0.38203125000000004, 1.0, 0.4860335195530726, 0.1655313351498638, 0.4644171779141105], 'peripheral artery disease': [0.5], 'Angelman syndrome': [0.5677233429394812, 0.9970760233918129, 0.09558823529411764], 'Facioscapulohumeral muscular dystrophy': [0.8092485549132948, 0.10755813953488375], 'epithelial ovarian cancer': [0.888888888888889, 0.9065420560747663, 0.9484240687679083, 0.8382066276803118, 0.9548104956268222], 'liver injury': [0.5], 'lung cancer': [0.4783498759305211, 0.2978868438991138, 0.6094377510040161, 0.3953373015873015, 0.3481583563016788], 'infertility': [0.27348066298342544, 0.1855345911949685], 'endometriosis': [0.17329545454545459, 0.9917355371900827, 0.5215384615384615], 'Heart Failure': [0.14244186046511625, 0.11490683229813664], 'Silver-Russell syndrome': [0.6536312849162011, 0.9814241486068112], 'Stroke': [0.5], 'Split Hand/Split Foot malformation disorder': [0.5], 'coronary artery disease': [0.8333333333333334, 0.940677966101695], "Parkinson's disease": [0.08912386706948638, 0.1818181818181818, 0.10217983651226159], 'acute myeloid leukemia': [0.7931034482758621, 0.8545706371191135, 0.9696048632218845, 0.20621468926553677, 0.9232954545454545]}
# ('AUC fold CV:', 0.669594891215425)
# # AUC for all diseases in CV:
# ('Burkitts lymphoma', 0.21549596845649477)
# ('esophageal adenocancer', 0.9170608548071348)
# ('pancreas cancer', 0.4192375174486777)
# ('nasopharyngeal cancer', 0.4245784278141162)
# ('Autism spectrum disorder', 0.5)
# ('panic disorder', 0.5)
# ('acute lymphoblastic leukemia', 0.27392707066989647)
# ('Intrauterine Growth Restriction', 0.5)
# ('B-cell neoplasms', 0.7903883724188773)
# ("Huntington's disease", 0.49582416643466487)
# ('B-cell lymphoma', 0.33921214725684007)
# ('dilated cardiomyopathy', 0.6652301713830211)
# ('rhabdomyosarcoma', 0.5)
# ('polymyositis', 0.5)
# ('Leishmania', 0.5)
# ('basal-like breast cancer', 0.5)
# ('promyelocytic leukemia', 0.5)
# ('pre-eclampsia', 0.5)
# ('atherosclerosis', 0.6656212984876304)
# ('multiple myeloma', 0.7392001553400539)
# ('gastric cancer', 0.403567085747659)
# ('Obesity', 0.3504911505286318)
# ('lymphoma', 0.5798960005620211)
# ('Wiedemann-Beckwith syndrome', 0.5)
# ('oesophageal adenocancer', 0.5)
# ('adolescent idiopathic scoliosis', 0.1375908641231881)
# ('DiGeorge syndrome', 0.5)
# ('postmenopausal osteoporosis', 0.5)
# ('chronic lymphocytic leukemia', 0.41367306521649816)
# ('autoimmune disease', 0.3051909730418772)
# ('myocardial infarction', 0.899512633063999)
# ('Wilms tumor', 0.9198736843549433)
# ('gestational trophoblastic diseases', 0.5)
# ('renal cell cancer', 0.3674988128123283)
# ('Crohn disease', 0.5)
# ('Duchenne muscular dystrophy', 0.4152275589134736)
# ('acute myocardial infarction', 0.15167822788800994)
# ('ovarian cancer', 0.6191349775164434)
# ('bone diseases', 0.5)
# ('schizophrenia', 0.8336653467920072)
# ('laryngeal cancer', 0.5)
# ("Hodgkin's lymphoma", 0.5)
# ('frontotemporal lobar degeneration', 0.5)
# ("Alzheimer's disease", 0.6066512420434098)
# ('myeloproliferative polycythaemia vera', 0.5)
# ('bipolar disorder', 0.9902961332327295)
# ('diabetic nephropathy', 0.5)
# ('basal cell cancer', 0.5)
# ('dyskeratosis congenita', 0.13354606168990454)
# ('psychiatric disease', 0.23077330077330077)
# ('cervical cancer', 0.6019237942880947)
# ('gallbladder cancer', 0.6553696716528761)
# ('pheochromocytoma', 0.9419061302681992)
# ('HCV', 0.5)
# ('Prader-Willi syndrome', 0.3249820730047217)
# ('acute megakaryoblastic leukemia', 0.9754572209375657)
# ('Neurodevelopmental syndromes associated with the SOX2 locus', 0.5)
# ('velocardiofacial syndrome', 0.5)
# ('neural tube defects', 0.5)
# ('cardiomyopathy', 0.5)
# ('endometrial cancer', 0.2708032214838824)
# ('type 2 diabetes', 0.4903392777858911)
# ('medulloblastoma', 0.5)
# ('adrenocortical cancers', 0.5)
# ('meningioma', 0.5)
# ('hepatocellular cancer', 0.3684597236586337)
# ('growth restriction', 0.5)
# ('Congenital hyperinsulinism', 0.5)
# ('testicular cancer', 0.2594704019754034)
# ('aortic aneurysm', 0.5)
# ('intracranial aneurism', 0.2998578483332101)
# ('cancer', 0.6197793710513076)
# ('malignant pleural mesothelioma', 0.4869870881605399)
# ('AIDS', 0.6677689422771532)
# ('trophoblastic tumor', 0.5)
# ('pseudohypoparathyroidism type Ib', 0.5)
# ('malignant melanoma', 0.707984366028148)
# ('neurofibromatosis type 1', 0.5)
# ('acute promyelocytic leukemia', 0.6020238095238095)
# ('gastric adenocancer', 0.33196277256233697)
# ('hepatoma', 0.5)
# ('biparental complete hydatidiform moles', 0.5)
# ('gastrointestinal stromal tumor', 0.5)
# ('chronic myeloid leukemia', 0.5)
# ('breast cancer', 0.5032154163954952)
# ('amyotrophic lateral sclerosis', 0.5)
# ('Diabetes', 0.6836720787249349)
# ('psoriasis', 0.6476189873844701)
# ('fibrosarcoma', 0.5)
# ('Glaucoma', 0.5)
# ('depression', 0.5)
# ('colon cancer', 0.4839481227097234)
# ('Pituitary adenoma', 0.714275063243548)
# ('oesophagus cancer', 0.5)
# ('chronic myeloproliferative disorders', 0.5)
# ('hematopoiesis', 0.5)
# ('small-cell lung cancer', 0.8259647371922821)
# ('Esophageal squamous cell cancer', 0.5073177741210623)
# ('myopia', 0.5)
# ('membranous nephropathy', 0.5)
# ('cleft lip', 0.5)
# ('squamous cell cancer', 0.5607670241445478)
# ('coronary disease', 0.5)
# ('papillary thyroid cancer', 0.28350305620325694)
# ('prostate cancer', 0.5193063872179924)
# ('Inclusion body myositis', 0.5)
# ('osteosarcoma', 0.4013671873862788)
# ('glioblastoma', 0.5737522582594756)
# ('hyperhomocysteinemia', 0.5)
# ('sarcoma', 0.5)
# ('tongue cancer', 0.9557700627298291)
# ('Beckwith-Wiedemann syndrome', 0.9068459333631559)
# ('periodontitis', 0.9742604497647726)
# ('lymphocytic leukemia', 0.5)
# ('aging', 0.35608308605341243)
# ('melanoma', 0.6772703021998142)
# ('hepatocelluar cancer', 0.6887902491615676)
# ('bladder cancer', 0.4913930455833463)
# ('drug abuse', 0.5)
# ('acute T-lymphocytic leukemia', 0.5)
# ('rheumatoid arthritis', 0.5)
# ('myelodysplastic syndrome', 0.24864584210378604)
# ('retinoblastoma', 0.5)
# ('multiple sclerosis', 0.5)
# ('coronary heart disease', 0.5)
# ('neuroblastoma', 0.49960265652340946)
# ('peripheral artery disease', 0.5)
# ('Angelman syndrome', 0.5534625338751372)
# ('Facioscapulohumeral muscular dystrophy', 0.45840334722408926)
# ('epithelial ovarian cancer', 0.9073744274077395)
# ('liver injury', 0.5)
# ('lung cancer', 0.4258340257445262)
# ('infertility', 0.22950762708919697)
# ('endometriosis', 0.5621898177579996)
# ('Heart Failure', 0.12867434638162645)
# ('Silver-Russell syndrome', 0.8175277167615062)
# ('Stroke', 0.5)
# ('Split Hand/Split Foot malformation disorder', 0.5)
# ('coronary artery disease', 0.8870056497175142)
# ("Parkinson's disease", 0.12437396179997658)
# ('acute myeloid leukemia', 0.7493578184855704)

# 10 fold
# ('AUC fold CV:', 0.6710498791710007)