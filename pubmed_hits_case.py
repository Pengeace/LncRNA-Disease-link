import pandas
import threading
import requests
from bs4 import BeautifulSoup


result_dir = './pubmed'

# base_url = 'https://www.ncbi.nlm.nih.gov/pubmed?term='

base_url = 'https://www.ncbi.nlm.nih.gov/pmc/?term='

diseases = ['Breast cancer']
lncRNAs = ["ACTA2-AS1","TINCR","CBR3-AS1","TUSC8","BDNF-AS","LINC00271","LINC00900","LINC00638","AC062029.1","ITGA9-AS1","RP11-66B24.4","C6orf3","LL0XNC01-116E7.2","LINC01004","MORF4L2-AS1","RP4-583P15.10","LINC00900","LINC00638","AC062029.1","ITGA9-AS1","RP11-66B24.4","C6orf3","LL0XNC01-116E7.2","LINC01004","MORF4L2-AS1","RP4-583P15.10","PSMD5-AS1"]


def get_pubmed_hits(disease):
    print('# disease: ' + disease)

    for lncRNA in lncRNAs:

        # search_target = str(lncRNA + ' AND ' + disease).replace(' ', '%20')
        search_target = str('(' + lncRNA + ')' + ' AND ' + '(' + disease + ')').replace(' ', '+')
        url = base_url+search_target
        page = requests.get(url)
        soup = BeautifulSoup(page.text, 'lxml')
        hit_count = int(soup.find(id='resultcount').attrs['value'])

        print("{}\t{}".format(lncRNA, hit_count))


for disease in diseases:

    thread = threading.Thread(target=get_pubmed_hits,args=(disease,))
    thread.start()


