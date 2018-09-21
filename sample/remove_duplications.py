import pandas as pd
import numpy as np

lncRNA_disease_path = '../data/lncRNA_disease.csv'
microRNA_disease_path = '../data/microRNA_disease.csv'
microRNA_lncRNA_path = '../data/microRNA_lncRNA.csv'

duplicated_lncRNAs = [['DANCR', 'DANCr'],['Xist', 'XIST'],['BACE1AS', 'BACE1-AS'],['RP11-672F9.1', 'rP11-672F9.1'],['HOTAIrM1', 'HOTAIRM1'],['TDRG1', 'TDrG1'],['DGCr5', 'DGCR5'],['CBR3-AS1', 'CBr3-AS1'],['Neat1', 'NEAT-1', 'NEAT1'],['Loc554202', 'LOC554202'],['LincRNA-p21', 'lincRNA-p21'],['UBE3A-ATS', 'Ube3a-ATS'],['TINCR', 'TINCr'],['PCAT1', 'PCAT-1'],['TRAF3IP2-AS1', 'TrAF3IP2-AS1'],['FENDRR', 'FENDrr'],['MIR155HG', 'MIr155HG'],['MIR17HG', 'MIr17HG'],['SOX2OT', 'SOX2-OT'],['Dnm3os', 'DNM3OS'],['lncRNA-LALR1', 'LncRNA-LALR1'],['CRNDE', 'CrNDE'],['member 1 opposite strand/antisense transcript 1 (KCNQ1OT1)', 'KCNQ1OT1'],['prostate cancer gene expression marker 1 (PCGEM1)', 'PCGEM1'],['HOX transcript antisense RNA (HOTAIR)', 'HOTAIR', 'HOTAIr'],['tumour suppressor candidate 7 (TUSC7)', 'TUSC7'],['cyclin-dependent kinase inhibitor 2B antisense RNA 1 (ANRIL)', 'ANRIL']]
duplicated_microRNAs = [['hsa-mir--196a-2', 'hsa-mir-196a-2']]
duplicated_diseases = [['Autism spectrum disorder', 'autism spectrum disorder'],['Esophageal squamous cell cancer', 'esophageal squamous cell cancer', 'Esophageal squamous cell- cancer'],['Obesity', 'obesity'],['Malignant mesothelioma', 'malignant mesothelioma'],['Fragile X syndrome', 'fragile X syndrome'],['Intrahepatic cholangiocarcinoma', 'intrahepatic cholangiocarcinoma'],['Pituitary adenoma', 'pituitary adenoma'],['small-cell lung cancer', 'small cell lung cancer'],['Heart Failure', 'heart failure'],['Facioscapulohumeral muscular dystrophy', 'facioscapulohumeral muscular dystrophy'],['B-cell lymphoma', 'b-cell lymphoma', 'B cell lymphoma'],['esophageal squamous cell cancer', 'Esophageal squamous cell- cancer'],['Hepatocellular carcinoma', 'hepatocellular carcinoma'],['Spinocerebellar ataxia type 7', 'spinocerebellar ataxia type 7']]


ld = pd.read_csv(lncRNA_disease_path)
md = pd.read_csv(microRNA_disease_path)
ml = pd.read_csv(microRNA_lncRNA_path)


microRNAs_bef = list(md['Name'])
lncRNAs_bef = list(ld['Name'])
diseases_bef = list(ld.columns[1:])
print(len(microRNAs_bef), len(lncRNAs_bef), len(diseases_bef))
# (677, 1831, 281)

def get_name_map(duplicated_pairs, name_map):
    for pair in duplicated_pairs:
        for x in pair:
            if x not in name_map:
                name_map[x] = pair[0]
    return name_map

lncRNA_map = get_name_map(duplicated_lncRNAs, {})
for x in lncRNAs_bef:
    if x not in lncRNA_map:
        lncRNA_map[x] = x
microRNA_map = get_name_map(duplicated_microRNAs, {})
for x in microRNAs_bef:
    if x not in microRNA_map:
        microRNA_map[x] = x
disease_map = get_name_map(duplicated_diseases, {})
for x in diseases_bef:
    if x not in disease_map:
        disease_map[x] = x

microRNAs_new = [x for x in microRNAs_bef if microRNA_map[x]==x]
lncRNAs_new = [x for x in lncRNAs_bef if lncRNA_map[x]==x]
diseases_new = [x for x in diseases_bef if disease_map[x]==x]
print(len(microRNAs_new), len(lncRNAs_new), len(diseases_new))
# (676, 1802, 266)


ld_new = pd.DataFrame(np.zeros([len(lncRNAs_new), len(diseases_new)]), index=lncRNAs_new, columns=diseases_new)
md_new = pd.DataFrame(np.zeros([len(microRNAs_new), len(diseases_new)]), index=microRNAs_new, columns=diseases_new)
ml_new = pd.DataFrame(np.zeros([len(microRNAs_new), len(lncRNAs_new)]), index=microRNAs_new, columns=lncRNAs_new)

ld = ld.set_index('Name')
md = md.set_index('Name')
ml = ml.set_index('Name')

for x in ld.index:
    for y in ld.columns:
        if ld.ix[x,y]==1:
            ld_new.ix[lncRNA_map[x], disease_map[y]] = 1

for x in md.index:
    for y in md.columns:
        if md.ix[x,y]==1:
            md_new.ix[microRNA_map[x], disease_map[y]] = 1

for x in ml.index:
    for y in ml.columns:
        if ml.ix[x,y]==1:
            ml_new.ix[microRNA_map[x], lncRNA_map[y]] = 1

lncRNA_disease_path = '../data/lncRNA_disease.csv'
microRNA_disease_path = '../data/microRNA_disease.csv'
microRNA_lncRNA_path = '../data/microRNA_lncRNA.csv'


ld_new.to_csv(lncRNA_disease_path)
md_new.to_csv(microRNA_disease_path)
ml_new.to_csv(microRNA_lncRNA_path)