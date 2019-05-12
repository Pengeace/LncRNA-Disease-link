import threading

import pandas
import requests
from bs4 import BeautifulSoup

result_dir = './PMC'
top_limit = 100
base_url = 'https://www.ncbi.nlm.nih.gov/pmc/?term='
# diseases = ['lung cancer', 'breast cancer', 'colorectal cancer', 'prostate cancer', 'ovarian cancer', "Alzheimer's disease", 'colon cancer','gastric cancer',"Parkinson's disease", 'pancreas cancer', 'testicular cancer', 'B-cell lymphoma', 'Stroke', "Huntington's disease", 'bipolar disorder', 'hepatocellular cancer']
# https://www.ncbi.nlm.nih.gov/pmc/?term=(Breast+cancer)+AND+(LINC01004)

diseases = ['Angelman syndrome', 'bipolar disorder', 'periodontitis', 'neurofibromatosis type 1', 'pheochromocytoma',
            'Beckwith-Wiedemann syndrome', 'meningioma', 'myocardial infarction', 'atherosclerosis', 'Wilms tumor',
            'schizophrenia', 'gastric adenocancer', 'epithelial ovarian cancer', 'acute myeloid leukemia',
            'Pituitary adenoma', 'hepatocelluar cancer', 'multiple myeloma', 'intracranial aneurism',
            'small-cell lung cancer', 'Diabetes', 'malignant pleural mesothelioma'] + \
           ['lung cancer', 'breast cancer',
            'colorectal cancer',
            'prostate cancer',
            'ovarian cancer',
            "Alzheimer's disease",
            'colon cancer', 'gastric cancer',
            "Parkinson's disease",
            'pancreas cancer',
            'testicular cancer',
            'B-cell lymphoma', 'Stroke',
            "Huntington's disease",
            'bipolar disorder',
            'hepatocellular cancer']

lncRNA_disease_link_TBSI = pandas.read_csv('lncRNA_disease_TBSI.csv').set_index('Unnamed: 0')


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
        print(lncRNA + disease + ' PMC-hits:' + str(hit_count))
    result_file = result_dir + './top ' + str(top_limit) + ' lncRNAs related to ' + disease + '.csv'
    disease_data.to_csv(result_file)

completes = ['top 100 lncRNAs related to acute myeloid leukemia.csv', "top 100 lncRNAs related to Alzheimer's disease.csv", 'top 100 lncRNAs related to Angelman syndrome.csv', 'top 100 lncRNAs related to atherosclerosis.csv', 'top 100 lncRNAs related to B-cell lymphoma.csv', 'top 100 lncRNAs related to Beckwith-Wiedemann syndrome.csv', 'top 100 lncRNAs related to bipolar disorder.csv', 'top 100 lncRNAs related to breast cancer.csv', 'top 100 lncRNAs related to colorectal cancer.csv', 'top 100 lncRNAs related to Diabetes.csv', 'top 100 lncRNAs related to epithelial ovarian cancer.csv', 'top 100 lncRNAs related to gastric adenocancer.csv', 'top 100 lncRNAs related to gastric cancer.csv', 'top 100 lncRNAs related to hepatocelluar cancer.csv', 'top 100 lncRNAs related to hepatocellular cancer.csv', "top 100 lncRNAs related to Huntington's disease.csv", 'top 100 lncRNAs related to intracranial aneurism.csv', 'top 100 lncRNAs related to lung cancer.csv', 'top 100 lncRNAs related to malignant pleural mesothelioma.csv', 'top 100 lncRNAs related to meningioma.csv', 'top 100 lncRNAs related to multiple myeloma.csv', 'top 100 lncRNAs related to myocardial infarction.csv', 'top 100 lncRNAs related to neurofibromatosis type 1.csv', 'top 100 lncRNAs related to pancreas cancer.csv', "top 100 lncRNAs related to Parkinson's disease.csv", 'top 100 lncRNAs related to periodontitis.csv', 'top 100 lncRNAs related to pheochromocytoma.csv', 'top 100 lncRNAs related to Pituitary adenoma.csv', 'top 100 lncRNAs related to schizophrenia.csv', 'top 100 lncRNAs related to small-cell lung cancer.csv', 'top 100 lncRNAs related to Stroke.csv', 'top 100 lncRNAs related to testicular cancer.csv', 'top 100 lncRNAs related to Wilms tumor.csv']


for disease in diseases:
    if disease in lncRNA_disease_link_TBSI.columns  and  ('top 100 lncRNAs related to '+disease+'.csv' not in completes):
        thread = threading.Thread(target=get_pubmed_hits, args=(disease,))
        thread.start()
