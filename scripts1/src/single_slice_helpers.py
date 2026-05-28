import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc

from typing import Dict, List, Union

from scipy.sparse import csr_matrix
import smpath.method as spm



##  custom for preparation to GSVA

def safe_load_processed_data(SAMPLE, subfolder, data_centralized_folder)->ad.AnnData:
    """safely load mzs intensities (anndata) with only good quality ions"""
    adata_mz = ad.read_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))
    adata_mz = adata_mz[:,
               adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
               "mz"].astype(str).tolist()]

    return adata_mz


def read_clean_annotations(SAMPLE, pathways_types, subfolder,
                           data_centralized_folder)-> Dict[str, pd.DataFrame]:
    """provide the dictionary of annotations dataframes,
     where the keys are the pathways_types"""
    annots_dict = {}
    for pw_type in pathways_types:
        tmpfi = os.path.join(data_centralized_folder,
                             subfolder,
                             f'{SAMPLE}-{pw_type}-CleanAnnotations.tsv')
        annots_dict[pw_type] = pd.read_csv(tmpfi,
                                          sep='\t', header=0, index_col=None)

    return annots_dict


def safe_load_data_and_clean_annotations(SAMPLE,
                                         pathways_types,
                                         subfolder, data_centralized_folder) ->(ad.AnnData, pd.DataFrame):
    """
    safely load mzs intensities (anndata) and annotations (dataframe):
     only good quality ions AND only clean annotations, for both objects
    """
    adata_mz = safe_load_processed_data(SAMPLE, subfolder, data_centralized_folder)
    annots_dict = read_clean_annotations(SAMPLE, pathways_types, subfolder, data_centralized_folder)

    adequate_mz = set()
    for pw_type in pathways_types:
        adequate_mz = adequate_mz.union(
            set(annots_dict[pw_type]['mz'].astype(str).unique())
        )
    adata_mz = adata_mz[:, list(adequate_mz)].copy()

    return adata_mz, annots_dict







