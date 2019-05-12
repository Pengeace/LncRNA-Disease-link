import threading

import pandas
import requests
from bs4 import BeautifulSoup

result_dir = '../results/PMC-validation/'
top_limit = 100
lncRNA_disease_prediction_path = '../data/lncRNA-disease-prediction.csv'
base_url = 'https://www.ncbi.nlm.nih.gov/pmc/?term='

# https://www.ncbi.nlm.nih.gov/pmc/?term=(Breast+cancer)+AND+(LINC01004)

# simple test
diseases = ["breast cancer", "Alzheimer's disease", "prostate cancer"]

# diseases = ['Angelman syndrome', 'bipolar disorder', 'periodontitis', 'neurofibromatosis type 1', 'pheochromocytoma',
#             'Beckwith-Wiedemann syndrome', 'meningioma', 'myocardial infarction', 'atherosclerosis', 'Wilms tumor',
#             'schizophrenia', 'gastric adenocancer', 'epithelial ovarian cancer', 'acute myeloid leukemia',
#             'Pituitary adenoma', 'hepatocelluar cancer', 'multiple myeloma', 'intracranial aneurism',
#             'small-cell lung cancer', 'Diabetes', 'malignant pleural mesothelioma'] + \
#             ['lung cancer', 'breast cancer',
#             'colorectal cancer',
#             'prostate cancer',
#             'ovarian cancer',
#             "Alzheimer's disease",
#             'colon cancer', 'gastric cancer',
#             "Parkinson's disease",
#             'pancreas cancer',
#             'testicular cancer',
#             'B-cell lymphoma', 'Stroke',
#             "Huntington's disease",
#             'bipolar disorder',
#             'hepatocellular cancer']

lncRNA_disease_link_TBSI = pandas.read_csv(lncRNA_disease_prediction_path).set_index('Unnamed: 0')


def get_pubmed_hits(disease):
    print('# disease: ' + disease)
    disease_data = lncRNA_disease_link_TBSI.ix[:, disease]
    disease_data = disease_data.sort_values(axis=0, ascending=False)
    top_lncRNAs = disease_data.index[0:top_limit]
    disease_data = disease_data.ix[top_lncRNAs,].to_frame()
    disease_data['PMC-hits'] = 0
    for lncRNA in top_lncRNAs:
        search_target = str('(' + lncRNA + ')' + ' AND ' + '(' + disease + ')').replace(' ', '+')
        url = base_url + search_target
        page = requests.get(url)
        soup = BeautifulSoup(page.text, 'lxml')
        hit_count = int(soup.find(id='resultcount').attrs['value'])

        disease_data.ix[lncRNA, 'PMC-hits'] = hit_count
        print("[" + lncRNA + ", " + disease + "] " + ' PMC-hits:' + str(hit_count))
    result_file = result_dir + disease + ' related top ' + str(top_limit) + ' lncRNAs' + '.csv'
    disease_data.rename(columns={disease: 'Score'}, inplace=True)
    # disease_data.to_csv(result_file, index=True, index_label="LncRNA")


for disease in diseases:
    if disease in lncRNA_disease_link_TBSI.columns:
        thread = threading.Thread(target=get_pubmed_hits, args=(disease,))
        thread.start()
