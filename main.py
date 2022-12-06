import os, argparse

import scipy
import pandas as pd
import numpy as np
import matlab.engine
import math
from scipy import io
from utils import dataloading_utils, method_utils, classifier_utils
from sklearn.preprocessing import OneHotEncoder
from scipy.stats import mode

from pycaret.classification import *
import rpy2.robjects as robjects


import anndata

import random

# random.seed(7)
random.seed(21)  #

if __name__ == '__main__':

    data_dir = "D:/scCPEP/codes/datas"


    ## parse arguments
    import argparse

    parser = argparse.ArgumentParser(description="data preprocessing")

    parser.add_argument('--DataPath', help="Data file path (.csv), cells-genes matrix with cell unique barcodes "
                                             "as row names and gene names as column names.",
                        default=None)

    parser.add_argument('--LabelsPath', help="Data label file path (.csv), cells-genes matrix with cell unique barcodes "
                                             "as row names and 'cell.type' as column names.",
                        default=None)
    parser.add_argument('--CV_RDataPath', help="Cross validation RData file path (.RData).",
                        default=None)
    parser.add_argument('--OutputDir', help="Output directory defining the path of the exported file.",
                        default=None)
    parser.add_argument('--sample_seed', help="Downsample seed in combined individual effect",
                        default=0, type=int)

    parser.add_argument('--feature_selection',
                        default=["DV", "DD", "chisq", "BI", "Ftest"])


    parser.add_argument('--rejection', help='whether apply(1) prediction rejection or not(0)',
                        default=0, type=float)

    parser.add_argument('--lognorm', default=0, type=int, help="whether apply(1) lognorm or not(0)")

    parser.add_argument('--iter', default=1, type=int)

    args = parser.parse_args()


    result_dir = "D:/scCPEP/codes/tmp"

    # read the Rdata file
    robjects.r['load'](args.CV_RDataPath)

    nfolds = np.array(robjects.r['n_folds'], dtype='int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype='bool')
    col = np.array(robjects.r['col_Index'], dtype='int')
    col = col - 1
    test_ind = np.array(robjects.r['Test_Idx'])
    train_ind = np.array(robjects.r['Train_Idx'])

    # read the data
    data = pd.read_csv(args.DataPath, index_col=0, sep=',')
    labels = pd.read_csv(args.LabelsPath, header=0, index_col=None, sep=',', usecols=col)

    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]

    for fold in range(np.squeeze(nfolds)):
    # for fold in (3, 5):
        test_ind_fold = np.array(test_ind[fold], dtype='int') - 1
        train_ind_fold = np.array(train_ind[fold], dtype='int') - 1

        df_train = data.iloc[train_ind_fold]
        df_train = df_train.T
        df_test = data.iloc[test_ind_fold]
        df_test = df_test.T
        obs_train = labels.iloc[train_ind_fold]
        obs_test = labels.iloc[test_ind_fold]

        obs_train.columns = ["cell.type"]
        obs_train.index = df_train.columns
        obs_test.columns = ["cell.type"]
        obs_test.index = df_test.columns

        var_train = pd.DataFrame(data=df_train.index, index=df_train.index)
        var_train.columns = ['gene_symbols']
        var_test = pd.DataFrame(data=df_test.index, index=df_test.index)
        var_test.columns = ['gene_symbols']

        train_adata = anndata.AnnData(X=df_train.T, obs=obs_train, var=var_train)
        train_adata.obs_names_make_unique(join="-")
        train_adata.var_names_make_unique(join="-")
        test_adata = anndata.AnnData(X=df_test.T, obs=obs_test, var=var_test)
        test_adata.obs_names_make_unique(join="-")
        test_adata.var_names_make_unique(join="-")

        
        train_adata, test_adata = dataloading_utils.preprocessing_data(train_adata, test_adata, result_dir, args=args)

        for it in range(args.iter):
            tmp_adata = train_adata.copy()
            ## to csv
            if scipy.sparse.issparse(tmp_adata.X) or \
                    isinstance(tmp_adata.X, pd.DataFrame):
                tmp_data = tmp_adata.X.toarray()
            else:
                tmp_data = tmp_adata.X

            ## write out original read count matrix
            tmp_df = pd.DataFrame(data=tmp_data, index=tmp_adata.obs_names, columns=tmp_adata.var_names).T
            tmp_df_path = result_dir + os.sep + "sub_tmp_counts.csv"
            tmp_df.to_csv(tmp_df_path)

            cell_annots = tmp_adata.obs["cell.type"].tolist()
            cell_annots_path = result_dir + os.sep + "sub_tmp_cell_annots.txt"
            with open(cell_annots_path, 'w') as f:
                for cell_annot in cell_annots:
                    f.write("%s\n" % cell_annot)
                    f.flush()
            f.close()
            del tmp_adata
            model_count = 0 # index of model

            celltype_cols = "cell.type"
            # OneHotEncoding the celltypes
            enc_train = OneHotEncoder(handle_unknown='ignore')
            #
            y_train = enc_train.fit_transform(train_adata.obs[[celltype_cols]]).toarray()
            y_test = test_adata.obs[[celltype_cols]]

            feature_selection = args.feature_selection
            i = 1
            for sel in feature_selection:

                ## each iter use one feature selection method
                train_sub, test_sub = method_utils.sub_feature_selection(sel,
                                                                         train_adata, test_adata,
                                                                         result_dir, tmp_df_path,
                                                                         cell_annots_path)
                if scipy.sparse.issparse(train_sub.X):
                    x_train = train_sub.X.toarray()
                else:
                    x_train = train_sub.X

                if scipy.sparse.issparse(test_sub.X):
                    x_test = test_sub.X.toarray()
                else:
                    x_test = test_sub.X

                trainXY = train_sub.to_df()
                trainXY['cell.type'] = y_train.argmax(1)
                testX = test_sub.to_df()

                if i == 1:
                    predictClass, testPred, model_count = classifier_utils.get_classifiers(args,
                                                                                           trainXY, testX,
                                                                                           data_dir, result_dir,
                                                                                           model_count)
                else:
                    pC, tP, model_count = classifier_utils.get_classifiers(args,
                                                                           trainXY, testX,
                                                                           data_dir, result_dir,
                                                                           model_count)
                    predictClass = np.append(predictClass, pC, axis=1)
                    testPred = np.append(testPred, tP, axis=1)
                i = i + 1
            trueClass = y_train.argmax(1)

            mat_path = "D:/scCPEP/codes" + os.sep + "ensemble.mat"

            scipy.io.savemat(mat_path, {'PredictClass': predictClass, 'TrueClass': trueClass})

            os.system("rm {}".format(cell_annots_path))  ## remove the temporaty cell annotations
            os.system("rm {}".format(tmp_df_path))  ## remove the temporaty counts

            print("=====here=====")

            # # PEP(MATLAB)
            eng = matlab.engine.start_matlab()
            # results = eng.pyTrain()
            results = eng.startTrain()
            ensembleResults = np.asarray(results[0], dtype=int)
            print("ensembleResults:", ensembleResults)
            ensemble = pd.DataFrame(ensembleResults)
            ensemble.to_csv(args.OutputDir + os.sep + "ensembleResult_" + str(fold) + ".csv", index=False) ## save ensemble resulr


            res_ind = np.where(ensembleResults == 1)
            tag = 0

            predResult = testPred[:, np.where(ensembleResults == 1)[0]]
            # # PEP result
            pred = mode(predResult.T)[0]



            ###### one result
            test_index = test_adata.obs.index
            p = pd.DataFrame(pred.T[:, 0], columns=[0], index=test_index)

            common_celltypes = set(train_adata.obs["cell.type"]).intersection(set(test_adata.obs["cell.type"]))
            test_cells = test_adata.obs.index
            

            prd = np.array(p.loc[test_cells])[:, 0]
            cellTypes_train = enc_train.categories_[0]
            add_unassigned = np.append(cellTypes_train, 'Unknown')
            prd_type = add_unassigned[prd.astype(np.int_)]
            
            trueLabs = pd.DataFrame(y_test.to_numpy())
            predLabs = pd.DataFrame(prd_type)
            
        trueLabs.to_csv(args.OutputDir + os.sep + "trueLabels_" + str(fold) + ".csv", index=False)
        predLabs.to_csv(args.OutputDir + os.sep + "predLabels_" + str(fold) + ".csv", index=False)
        ########## one result

