# scCPEP
scCPEP method aims at the cell-type annotation problem. This method first
takes known cell types as a reference and then identifies
the cell types of reference using ensemble models created
by integrating multiple feature selection approaches and
classifiers. Then, an optimal subset of the ensemble models
is extracted using Pareto Ensemble Pruning (Qian, Yu and
Zhou, 2015) technique. This subset of the models is then
used to identify the specific cell types in the query dataset,
and the prediction results of the subset are aggregated to
provide a final prediction. 
