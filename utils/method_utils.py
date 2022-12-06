import os, sys
import numpy as np
import pandas as pd
import scipy
import scanpy as sc
import anndata

from sklearn import metrics





## ---- some functions for processing data
def process_adata(adata, min_genes=10, min_cells=10, lognorm=True, celltype_label="cell.type"):
    '''Procedures for filtering single-cell data
       1. Filter low-quality cells and genes;
       2. Filter nonsense genes;
       3. Normalize and log-transform the data;
       4. Change all gene names into UPPER;
       5. Remove cells with no labels;
    '''
    adata.var_names = [i.upper() for i in list(adata.var_names)]#avod some genes having lower letter

    ## make names unique
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # pre filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)

    # pre_filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # prefilter_specialgene: MT and ERCC  -> from ItClust package
    Gene1Pattern="ERCC"
    Gene2Pattern="MT-"
    id_tmp1 = np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names], dtype=bool)
    id_tmp2 = np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names], dtype=bool)
    id_tmp = np.logical_and(id_tmp1, id_tmp2)
    adata._inplace_subset_var(id_tmp)

    ## handel exception when there are not enough cells or genes after filtering
    if adata.shape[0] < 3 or adata.shape[1] < 3:
        return None

    #4 normalization,var.genes,log1p
    if lognorm == 1:
        sc.pp.normalize_per_cell(adata)  ## total count equal to the median of the counts_per_cell
        sc.pp.log1p(adata)

    ## cells with celltypes
    cells = adata.obs.dropna(subset=[celltype_label]).index.tolist()
    adata = adata[cells]
    return adata


def preprocess_train_test(train_adata, test_adata, result_dir,
        gene_no=1000, lognorm=True, select_on="test", select_method="Seurat",
        min_genes=10, min_cells=10, celltype_label="cell.type"):

    ## filter cells/genes
    train_adata = process_adata(train_adata, min_genes, min_cells, lognorm)
    test_adata = process_adata(test_adata, min_genes, min_cells, lognorm)

    ## handle None exception
    if train_adata is None or test_adata is None:
        return None, None
    return train_adata, test_adata

    


DV_PATH = "Rimpl/doDV.R"
limma_PATH = "Rimpl/doLimma.R"
DD_PATH = "Rimpl/doDD.R"
chisq_PATH = "Rimpl/doChisSquared.R"
BI_PATH = "Rimpl/doBI.R"
FTEST_RSCRIPT_PATH = "Rimpl/Ftest_selection.R"
def sub_feature_selection(method, train_adata, test_adata,
                          result_dir, tmp_df_path, cell_annots_path,
                          celltype_label="cell.type"):
    topN = 50
    pSig = 0.001
    gene_no = 2000

    # if "limma" == method:
        # os.system("Rscript --vanilla " + limma_PATH + " " + tmp_df_path + " " +
        #           cell_annots_path)
    if "Ftest" == method:
        os.system("Rscript --vanilla " + FEAST_FTEST_RSCRIPT_PATH + " " + tmp_df_path + " " +
                  cell_annots_path + " " + str(gene_no))
    if "DV" == method:
        os.system("Rscript --vanilla " + DV_PATH + " " + tmp_df_path + " " +
                  cell_annots_path)
    if "DD" == method:
        os.system("Rscript --vanilla " + DD_PATH + " " + tmp_df_path + " " +
                  cell_annots_path)

    if "chisq" == method:
        os.system("Rscript --vanilla " + chisq_PATH + " " + tmp_df_path + " " +
                  cell_annots_path)
    if "BI" == method:
        os.system("Rscript --vanilla " + BI_PATH + " " + tmp_df_path + " " +
                  cell_annots_path)
                  
                  
    ## handle None exception
    if train_adata is None or test_adata is None:
        return None, None

    ## read features selected by each method
    sub_feature_file = result_dir + os.sep + "sub_features.txt"
    with open(sub_feature_file) as f:
        features = f.read().splitlines()
    sub_features = [x.upper() for x in features]  ## upper case

    sub_genes = set(sub_features).intersection(set(train_adata.var_names.tolist()))
    train_adata = train_adata[:, list(sub_genes)]

    features = set(train_adata.var_names.tolist()).intersection(set(test_adata.var_names.tolist()))
    features = list(features)
    features.sort()  ## for reproducibility
    print("Number of", method, "features:", len(features))

    ## write features into file
    with open(result_dir + os.sep + "sub_features.txt", 'w') as f:
        for feature in features:
            f.write("%s\n" % feature)
            f.flush()
    f.close()

    ## order common genes in anndata
    train_adata = train_adata[:, features]
    test_adata = test_adata[:, features]
    return train_adata, test_adata



