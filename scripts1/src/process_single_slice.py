from typing import  Tuple
import logging
import warnings
warnings.simplefilter('ignore',category=FutureWarning)
import numpy as np
from smpath.data import MetabolomeAnnotations
from smpath.processing.measures_methods import compute_moran_i
from smpath.helpers import setup_logger  # to add extra info to the smpath log
import smpath.data as spd
from sklearn.cluster import KMeans



from .single_slice_helpers import *

logger = logging.getLogger(__name__)
setup_logger(logger)


np.random.seed(42)  # important to make the computations be reproducible


pathways_files_location_info = {
        "KEGG" : {'file': "KEGG_hsa_pathways_compounds_R110.gmt",
                      'subfolder': "metabolites-mappings",
                      'id_type': 'KEGG'},
        "WIKI" : {'file': "WIKIPATHWAYS_human_multi_id_redundant.mmt",
                          'subfolder': "metabolites-mappings",
                          'id_type': 'moleculeIds'}}


def join_clean_annotations(pathways_types: List[str],
                    SAMPLE: str, subfolder: str,
                   data_centralized_folder: str)-> pd.DataFrame:
    """
    concatenates annotations for each element in pathways_types
     (but not FAMILIES, as the annotations in this case are the same that served to 'WIKI' type of pathway)
     in other words, concatenates KEGG and no-KEGG (that we coined 'WIKI') annotations
    """
    pathways_types = list(set(pathways_types) - {"FAMILIES"})
    annots_dict = read_clean_annotations(SAMPLE, pathways_types, subfolder, data_centralized_folder)
    dfs_to_gather = list()
    for k in annots_dict.keys():
        if k == 'KEGG':
            df = annots_dict[k][["KEGG", "moleculeNames"]]
            df.columns = ["moleculeIds", "moleculeNames"]
        elif k == 'WIKI':
            df = annots_dict[k][["moleculeIds", "moleculeNames"]]
        dfs_to_gather.append(df)

    gathered_annotations = pd.concat(dfs_to_gather, axis=0)
    return gathered_annotations


def compute_compounds_df_and_save(pathways_types, min_pathway_size, SAMPLE,
                           subfolder, data_centralized_folder):
    for pw_type in pathways_types:
        adata_mz, clean_annots_dict = safe_load_data_and_clean_annotations(
            SAMPLE, pathways_types, subfolder, data_centralized_folder
        )
        annots_df = clean_annots_dict[pw_type]
        pathways_obj = spd.MetabolicPathways(
            data_dir=data_centralized_folder,
            subfolder=pathways_files_location_info[pw_type]['subfolder'],
            file=pathways_files_location_info[pw_type]['file']
        )
        pathways_obj.molecule_id_type = pathways_files_location_info[pw_type]['id_type']

        pathways_obj.load()
        pathways_obj, annots_df = spm.SyncPathwaysAndAnnots.run(
            pathways_obj.model_copy(),
            annots_df,
            usable_id=pathways_files_location_info[pw_type]['id_type'],
            min_set_size=min_pathway_size,
        mode="mz_aware"  #
        )

        metabolites_computed_df, coords = spm.ComputeAssignmentsWeights.run(
            adata_mz, annots_df,
            usable_id=pathways_files_location_info[pw_type]['id_type'])

        metabolites_computed_df.to_csv(
            os.path.join(data_centralized_folder,
                                        subfolder,
                         f'{SAMPLE}-cmps-mx-{pw_type}.tsv'),
            sep='\t', index=True, header=True)

    # coords outside loop as any one has same coords
    coords_df = pd.DataFrame(data=coords,
                             index=metabolites_computed_df.index.tolist(),
                             columns=['x', 'y'])
    coords_df.to_csv(os.path.join(data_centralized_folder,
                                    subfolder, f'{SAMPLE}-coords.tsv'),
        sep='\t', index=True, header=True)


# -----------------------
# ######### end #########
# not used:
#
# def table_contracted(df):
#     df = df.drop(columns=['moleculeIds'])
#     transformed_df = df.groupby(['mz', 'cluster', 'pathway_name',
#                                  'pathway_id', 'theoretical_direction'], as_index=False).agg({
#         'moleculeNames': lambda x: '; '.join(list(set(x)))
#     })
#     return transformed_df
