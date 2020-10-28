import pandas as pd
import os


"""
Initial parameters set through the configuration file.
"""
filename_config = 'config.ini'
pathway_input = None
gene_input = None
deep_input = None
num_cores_input = None
mode_input = None


"""
Parameters set during execution.
"""
gene_input_hsa = None
gene_input_url = None


"""
display.max_colwidth sets the maximum width of columns.
If the value is "None", the max length is disabled.
"""
pd.set_option('display.max_colwidth', None)


"""
At this point, the dataframe with the specified columns will be created."""
COLS_DF = ['deep', 'name_father', 'hsa_father', 'name_son', 'hsa_son', 'url_kegg_son', 'relation',
           'type_rel', 'pathway_of_origin', 'fullpath', 'occurrences']
DF_TREE = pd.DataFrame(columns=COLS_DF)

COLS_DF_2 = ['hsa', 'gene']
DF_GENE_HSA = pd.DataFrame(columns=COLS_DF_2)


"""
At this point, the dictionary will be created and used to generate the json file.
"""
json_dict = {}
