- *create_index.py* gives the index number to each lncRNA, microRNA and disease in the heterogeneous tripartite network.
- *deepwalk_training.py* generates the feature vector for each biomedical node in the tripartite network.
- *rule_based_prediction.py* calculates the association score for each lncRNA-disease pair.
- *PMC_hits.py* and *Pubmed_hits.py* are for prediction result validation on [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/) 
and [PubMed Center](https://www.ncbi.nlm.nih.gov/pmc/)(PMC).
- The rest of *deep_walk_training_CV.py*, *prediction_CV.py* and *parallel_prediction_CV.py* are for 5-fold cross validation (CV). 
