from .process_single_slice import *
import warnings
warnings.simplefilter('ignore',category=FutureWarning)



def annotation_types_to_file_names_dict(SAMPLE, annotation_types_list):
    """fills dictionary of annotation files with standardized file names"""
    out_d = {}
    for annotation_type in annotation_types_list:
        out_d[annotation_type] = f"{SAMPLE}-{annotation_type}.csv"

    return out_d



def process_tissues(input_file:str, data_dir:str
                     )-> None:
    """
    - parses the input dataframe having the specifics of the input datasets
        for each row representing a dataset.
    - computes the matrix of compounds (compounds as columns, spots as rows)
        with the probabilistic approach as in SPACEPATH
    - saves the matrix for comparison with other methodologies: GSVA
    """
    data_centralized_folder = data_dir

    datasets_parameters_df = pd.read_csv(input_file,
        header=0, index_col=None, sep='\t')

    for i, row in datasets_parameters_df.iterrows():

        SAMPLE = row["dataset_name"]
        annotations_to_use = row["annotations_to_use"].split(" ")

        pathways_types = row["pathways_types"].split(" ")
        min_pathway_size = row["min_pathway_size"]


        logger.info('\nRunning sample: %s', SAMPLE)

        subfolder = SAMPLE

        compute_compounds_df_and_save(pathways_types, min_pathway_size,
                                      SAMPLE,
                                      subfolder, data_centralized_folder)

        logger.info("processed sample %s \n\n", SAMPLE)
