import pandas as pd
import anndata
import os, sys
from utils import method_utils


def load_data(datapath, datalabel=None):
    df = pd.read_csv(datapath, index_col=0)
    df = df.T
    obs = pd.read_csv(datalabel, index_col=0)

    ## build obs/var dataframe
    obs.columns = ["cell.type"]
    obs.index = df.columns
    var = pd.DataFrame(data=df.index, index=df.index)
    var.columns = ['gene_symbols']

    adata = anndata.AnnData(X=df.T, obs=obs, var=var)
    adata.obs_names_make_unique(join="-")
    adata.var_names_make_unique(join="-")

    return adata




def preprocessing_data(train_adata, test_adata, result_dir, args=None, scale=True, plot=True, save_raw=False,
                       save_data=True):

    if args is None:
        sys.exit("Error: Please check your argument parser object!")

    if save_raw:
        train_adata.layers["counts"] = train_adata.X.copy()
        test_adata.layers["counts"] = test_adata.X.copy()

    # preprocessing
    train_adata, test_adata = method_utils.preprocess_train_test(train_adata, test_adata, result_dir)

    return train_adata, test_adata

