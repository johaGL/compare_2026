import os
import anndata as ad
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr
from scipy.spatial import distance
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc

from smpath.helpers import generate_coordinates_arr
import smpath.processing.measures_methods as spm
from scipy.stats import zscore
import squidpy as sq

palette_box = {'kpca' : "#66c2a5", 'gsva': "skyblue"}
palette_swarm = {'kpca': "green", 'gsva': "cadetblue" }


def open_scores_objects_full(out_parent_dir, dir_adata, SAMPLE):
    """
    opens both gsva and kpca matrices (dataframes)
    the columns' names are the set identifiers, the rows' names the spots.
    returns both separated matrices with the union of their column names (some columns are NaN)
    """
    matrix_gsva_kegg = pd.read_csv(os.path.join(
        out_parent_dir, SAMPLE, f"{SAMPLE}-cmps-mx-KEGG-KEGG_gsva.tsv"
    ), sep='\t', index_col=0)
    adata_score_kegg = ad.read_h5ad(os.path.join(dir_adata, SAMPLE,
                                                 f"{SAMPLE}-scores-KEGG.h5ad"))
    kpca_kegg = adata_score_kegg.to_df()
    print(
        f"initial sets numbers: kpca -> {kpca_kegg.shape[1]}, gsva -> {matrix_gsva_kegg.shape[1]}")

    ids_lack_in_a = set(matrix_gsva_kegg.columns.tolist()) - set(kpca_kegg.columns.tolist())
    ids_lack_in_b = set(kpca_kegg.columns.tolist()) - set(matrix_gsva_kegg.columns.tolist())

    ok_spots = list(set(kpca_kegg.index.tolist()
                        ).intersection(set(matrix_gsva_kegg.index.tolist())))

    kpca_kegg = kpca_kegg.loc[ok_spots, :]
    for col in ids_lack_in_a:
        kpca_kegg[col] = np.nan
    matrix_gsva_kegg = matrix_gsva_kegg.loc[ok_spots, :]
    for col in ids_lack_in_b:
        matrix_gsva_kegg[col] = np.nan

    return kpca_kegg, matrix_gsva_kegg


def subset_intersection_scores_objects(kpca_kegg:pd.DataFrame, matrix_gsva_kegg:pd.DataFrame):
    """
    opens both matrices (dataframes)
    drops NaN columns
    returns both separately, with the intersection of the columns' names
    """
    kpca_kegg = kpca_kegg.dropna(axis=1, how="all")
    matrix_gsva_kegg = matrix_gsva_kegg.dropna(axis=1, how="all")

    ok_ids = list(set(kpca_kegg.columns.tolist()
                      ).intersection(set(matrix_gsva_kegg.columns.tolist())))

    kpca_kegg = kpca_kegg.loc[:, ok_ids]
    matrix_gsva_kegg = matrix_gsva_kegg.loc[:, ok_ids]

    return kpca_kegg, matrix_gsva_kegg


# def open_scores_objects(out_parent_dir, dir_adata, SAMPLE):
#     matrix_gsva_kegg = pd.read_csv(os.path.join(
#         out_parent_dir, SAMPLE, f"{SAMPLE}-cmps-mx-KEGG-KEGG_gsva.tsv"
#     ), sep='\t', index_col=0)
#     adata_score_kegg = ad.read_h5ad(os.path.join(dir_adata, SAMPLE,
#                                                  f"{SAMPLE}-scores-KEGG.h5ad"))
#     kpca_kegg = adata_score_kegg.to_df()
#
#     print(
#         f"initial sets numbers: kpca -> {kpca_kegg.shape[1]}, gsva -> {matrix_gsva_kegg.shape[1]}")
#
#     print("  ", kpca_kegg.columns)
#     print("  ", matrix_gsva_kegg.columns)
#
#     ok_ids = list(set(kpca_kegg.columns.tolist()
#                       ).intersection(set(matrix_gsva_kegg.columns.tolist())))
#
#     ok_spots = list(set(kpca_kegg.index.tolist()
#                         ).intersection(set(matrix_gsva_kegg.index.tolist())))
#
#     kpca_kegg = kpca_kegg.loc[ok_spots, ok_ids]
#     matrix_gsva_kegg = matrix_gsva_kegg.loc[ok_spots, ok_ids]
#
#     return kpca_kegg, matrix_gsva_kegg


def return_pearsons_all_feats(kpca_kegg,
                              matrix_gsva_kegg) -> pd.DataFrame:
    kpca_kegg, matrix_gsva_kegg = subset_intersection_scores_objects(
        kpca_kegg, matrix_gsva_kegg
    )
    correls_res = list()

    assert np.all(
        kpca_kegg.columns == matrix_gsva_kegg.columns), "Error, columns do not match, aborting!"
    assert np.all(
        kpca_kegg.index == matrix_gsva_kegg.index), "Error, rows do not match, aborting!"

    for u in kpca_kegg.columns:
        res = pearsonr(matrix_gsva_kegg[u].to_numpy(),
                       kpca_kegg[u].to_numpy())
        correls_res.append(res.statistic)

    return pd.DataFrame({"pathway_id": kpca_kegg.columns, "r": correls_res})


def return_variances_of_scores(kpca_kegg: pd.DataFrame,
                               matrix_gsva_kegg: pd.DataFrame,
                               return_spotwise=False):
    assert np.all(
        kpca_kegg.columns == matrix_gsva_kegg.columns), "Error, columns do not match, aborting!"
    assert np.all(
        kpca_kegg.index == matrix_gsva_kegg.index), "Error, rows do not match, aborting!"

    axis_int = 0  # by default, across pathways (faster computation)
    if return_spotwise:
        axis_int = 1   # across spots (very slow when too many spots)

    flattened_kpca = kpca_kegg.to_numpy().var(axis=axis_int)  # .flatten()
    flattened_gsva = matrix_gsva_kegg.to_numpy().var(
        axis=axis_int)  # .flatten()

    types_list = list(np.repeat("kpca", len(flattened_kpca)
                                )) + list(
        np.repeat("gsva", len(flattened_gsva)))

    df = pd.DataFrame(
        {"avg_scores": list(flattened_kpca) + list(flattened_gsva),
         "type": types_list})

    return df


def compute_ordering_tissues(huge: pd.DataFrame, which_values: str):
    """
    huge is a dataframe with columns: pathway_id, <which_values>, tissue
    ranks the tissue names by the median of the values of the
        column <which_values>     across pathway_id s
    """
    assert which_values in ['r', 'abs(r)', 'value'], "Error, column not recognized"
    choo = pd.pivot(huge, columns="pathway_id", index="tissue",
                    values=which_values)
    choo = choo.assign(median_score_per_tissue=choo.median(axis=1).to_numpy())
    choo = choo.assign(tissue_name=choo.index.tolist())
    choo = choo[['median_score_per_tissue', 'tissue_name']]
    choo = choo.sort_values(by=['median_score_per_tissue'], ascending=True)

    return choo.index.tolist()


def do_binary_mask(v_array, centile_value=0.7):
    assert centile_value <= 1, "Error, centile_value must be between 0 and 1"
    try:
        value = np.percentile(v_array, q=centile_value * 100)
        print(v_array.max(), v_array.min())
        print("centile; ", value)
        mask = np.zeros(shape=v_array.shape, dtype=np.uint32)
        mask[mask >= value] = 1
    except Exception as e:
        print(e)

    return mask


def dice_coefficient(mask1, mask2):
    mask1 = np.asarray(mask1, dtype=bool)
    mask2 = np.asarray(mask2, dtype=bool)

    if mask1.shape != mask2.shape:
        raise ValueError("Masks must have the same shape")

    intersection = np.logical_and(mask1, mask2).sum()

    return 2.0 * intersection / (mask1.sum() + mask2.sum())


def plot_this_pathway(pathway_id, pdseries, modality: str,
                      ax):
    coords_strings_l = pdseries.index.tolist()
    if ('x' in coords_strings_l[0]) and ('y' in coords_strings_l[0]):
        coords_arr = generate_coordinates_arr(coords_strings_l)

    df = pd.DataFrame(data=coords_arr, columns=['x', 'y'])
    df = df.assign(values=pdseries.to_numpy())
    sns.scatterplot(data=df, x='x', y='y', hue='values',
                    palette='Spectral_r', legend=True,
                    s=5,
                    ax=ax)
    ax.set_title(f'{pathway_id} : {modality}')

    return ax


def wrap_comp_moran(df:pd.DataFrame, coords_df:pd.DataFrame,
                    set_type:str, tissue:str):
    # shift values to be all positive so avoids errors
    min_val_all = df.min(skipna=True).min(skipna=True)

    cols = df.columns.tolist()
    vars = df.index.tolist()

    df = df + abs(min_val_all) + 1e-20
    df = pd.DataFrame(df, index=vars, columns=cols)
    coords_df['obs'] = coords_df.index.tolist()
    coords_df.index = coords_df.index.astype(str).tolist()

    adata = ad.AnnData(
        X=csr_matrix(df.values),
        var=pd.DataFrame({'pathway_id': np.array(df.columns.tolist(), dtype=str)},
                         index=df.columns.tolist(), dtype=str),
        obs=coords_df,
    )
    adata.obsm['spatial'] = np.array(coords_df)
    adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 1,
                         'tissue_lowres_scalef': 1}}}
    try:
        sq.pl.spatial_scatter(adata, color=adata.var.index.tolist() ,img=False,
                              save=os.path.join("fofo", f'{tissue}--{set_type}.png'), dpi=200)

        print("yesssss!!!!!")

    except Exception as e:
        print("$$$$$$$$$$$$$$$$$===============", e)

    # adata = spm.compute_moran_i(adata, n_perms=2) # fast, no need pvalues here
    #
    # moran_df = adata.uns['moranI'].copy()
    # moran_df['pathway_id'] = moran_df.index.tolist()
    # moran_df['type'] = set_type
    # moran_df['tissue'] = tissue

    #return moran_df[['pathway_id', 'I', 'type', 'tissue']]
    return 0



if __name__ == '__main__':

    out_parent_dir = '../../apollo-gsva/compar_out_2026'  # TODO: modify in server: compar_out_2026
    dir_adata = "../../apollo-data" # TODO: modify in server: data

    tisues_df = pd.read_csv("../the_data_list.tsv", sep="\t",
                              index_col=None, header=0)
    tissues_list = tisues_df['dataset_name'].tolist()  #  [-4:]
    print(tissues_list)

    correls_dfs = list()
    scores_each_dfs = list()
    morans_dfs = list()

    for SAMPLE in tissues_list:
        try:
            adata_score_kegg = ad.read_h5ad(os.path.join(dir_adata, SAMPLE,
                                                         f"{SAMPLE}-scores-KEGG.h5ad"))

            coords_df = pd.read_csv(os.path.join(out_parent_dir,
                                                 SAMPLE, f'{SAMPLE}-coords.tsv'),
                                    sep='\t', index_col=0, header=0)
        except Exception as e:
            print(e, "failed opening")

      

    # -- now plot the moran I scores of the features for each gsva and kpca
    moran_comb = pd.concat(morans_dfs, ignore_index=True)

    moran_comb.to_csv("noentiendo.csv", index=False)

    # moran_comb = moran_comb.assign(
    #     tissue=pd.Categorical(moran_comb['tissue'].tolist(),
    #                           categories=tissues_order_by_correl))

    sns.boxplot(data=moran_comb, x="tissue", y="I", hue="type",
                palette=palette_box,
                legend=True,
                flierprops={"marker": "o"}, zorder=0)
    sns.stripplot(data=moran_comb, x="tissue", y="I", hue="type",
                  size=2.2, dodge=True,
                  palette=palette_swarm,
                  zorder=1,
                  legend=False)
    plt.xticks(rotation=90)
    plt.show()
    plt.close()


