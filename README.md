# LncRNA-Disease association prediction

Predicting lncRNA-disease associations using network topological similarity based on deep mining heterogeneous networks.

The programs are in Python 2.7. In [sample](https://github.com/Pengeace/LncRNA-Disease-link/tree/master/sample) 
and [data](https://github.com/Pengeace/LncRNA-Disease-link/tree/master/data) directories, the function of each script and data file is briefly described.

## Input

A heterogeneous tripartite network consists of three types of nodes, lncRNAs, microRNAs and diseases, 
were constructed from three kinds of bipartite networks.

- lncRNA-disease association network
- lncRNA-microRNA association network
- microRNA-disease associations network

## Method

Firstly, [DeepWalk](https://github.com/phanein/deepwalk) model was trained on the tripartite network to
generate the feature representation for each biomedical node in the network.

Then the similarity between every two lncRNAs was calculated based on the cosine distance of their feature vectors.

Finally, the association score of each lncRNA-disease pair was calculated by rule based method.

## Validation

The relation of predicted top lncRNA-disease pairs were validated by text mining in [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/) 
and [PubMed Center](https://www.ncbi.nlm.nih.gov/pmc/)(PMC). The whole validation results are in 
[result](https://github.com/Pengeace/LncRNA-Disease-link/tree/master/results) directory.

Below Table 1 listed the predicted Alzheimer's disease related top 5 lncRNAs and the corresponding
PMC hit counts. 

LncRNA | Score | PMC hit counts
----|----------|----------
51A | 0.02956689308352432 | 107
GDNFOS | 0.028450111719631625 | 24
HAR1A | 0.022572780455485025 | 11
HAR1B | 0.022258334177991525 | 3
HTTAS_v1 | 0.018553183446962912 | 8

For example, the below query found 107 records of Alzheimer's disease and 51A in PMC on March 12, 2019.

    https://www.ncbi.nlm.nih.gov/pmc/?term=Alzheimer's+disease+51A

---

[1] Perozzi, Bryan, Rami Al-Rfou, and Steven Skiena. "Deepwalk: Online learning of social representations." Proceedings of the 20th ACM SIGKDD international conference on Knowledge discovery and data mining. ACM, 2014.
