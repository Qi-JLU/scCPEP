import time, os
import tensorflow as tf
# import tensorflow.compat.v1 as tf
import numpy as np
import pandas as pd
# import h2o
from numpy.random import seed

#Set seeds
RANDOM_SEED=0
seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
tf.compat.v1.set_random_seed(RANDOM_SEED)

## every time reset all graphs and run something new
tf.compat.v1.reset_default_graph()
#tf.reset_default_graph()




## pycaret
import gc
from pycaret.classification import *
def get_classifiers(args,
                    train, testX,
                    data_dir, result_dir, model_count):
    '''Run methods
    '''

    from sklearn.naive_bayes import MultinomialNB
    from sklearn.svm import SVC, LinearSVC
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.naive_bayes import GaussianNB
    from sklearn.neural_network import MLPClassifier

    

    rejection = args.rejection

    model_num = 10
    
    threshold = rejection

    
    clf = setup(data=train, target='cell.type', fold=5,
                preprocess = False, silent = True)


    best = compare_models(n_select=model_num, include=[MultinomialNB(alpha=0.01),
                                                       CalibratedClassifierCV(LinearSVC()),
                                                       'mlp',
                                                       'knn',
                                                       'rbfsvm',
                                                       'rf'
                                                       ])


    for i in range(min(model_num, len(best))):
        model = best[i]
        # final_model = finalize_model(model)
        pred = predict_model(model, data=train)
        test_pred = predict_model(model, data=testX)
        if 0 == i:
            predictClass = pred[['Label']]
            testPred = test_pred[['Label']]

            ############### unassigned ###############c
            num_celltype = len(np.unique(train[['cell.type']].to_numpy()))
            testScore = test_pred[['Score']].to_numpy()
            unlabeled = np.where(testScore < threshold)
            # testPred[unlabeled[0]] = 'unassigned'
            testPred.iloc[unlabeled[0].tolist()] = num_celltype

            ## add 'unassigned' to predClass
            predScore = pred[['Score']].to_numpy()
            unlabeled_pred = np.where(predScore < threshold)
            predictClass.iloc[unlabeled_pred[0].tolist()] = num_celltype
            ##########################################
            # trueClass = pred[['cell.type']]
        else:
            # predictClass = np.append(predictClass, pred[['Label']], axis=1)
            # testPred = np.append(testPred, test_pred[['Label']], axis=1)

            pC = pred[['Label']]
            tp = test_pred[['Label']]

            ############## unassigned ###############
            num_celltype = len(np.unique(train[['cell.type']].to_numpy()))
            testScore = test_pred[['Score']]
            unlabeled = np.where(testScore < threshold)
            # testPred[unlabeled] = 'unassigned'
            tp.iloc[unlabeled[0].tolist()] = num_celltype
            testPred = np.append(testPred, tp, axis=1)

            ## add 'unassigned' to predClass
            predScore = pred[['Score']].to_numpy()
            unlabeled_pred = np.where(predScore < threshold)
            pC.iloc[unlabeled_pred[0].tolist()] = num_celltype
            predictClass = np.append(predictClass, pC, axis=1)
        #     ##########################################


    gc.collect()

    return predictClass, testPred, model_count

