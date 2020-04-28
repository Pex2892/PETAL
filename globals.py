import pandas as pd


"""
Initial parameters set through the configuration file.
"""
filename_config = 'config.json'
logger = None
pathway_input = None
gene_input = None
hop_input = None
num_cores_input = None


"""
Parameters set during execution.
"""
gene_input_hsa = None
gene_input_url = None


"""
display.max_colwidth sets the maximum width of columns.
If the value is -1, the max length is disabled.
"""
pd.set_option('display.max_colwidth', -1)


"""
At this point, the dataframe with the specified columns will be created."""
COLS_DF = ['hop', 'name_start', 'hsa_start', 'name_end', 'hsa_end', 'url_gene_end', 'relation',
           'type_rel', 'pathway_origin', 'path', 'occurrences']
DF_TREE = pd.DataFrame(columns=COLS_DF)


"""
At this point, the dictionary will be created and used to generate the json file.
"""
json_dict = {}


"""
This color dictionary is used to colorize prints during execution.
"""
COLORS = {
    "pink": "\033[95m",
    "blue": "\033[94m",
    "green": "\033[92m",
    "yellow": "\033[93m",
    "red": "\033[91m",
    "end_line": "\033[0m",
    "bold": "\033[1m",
    "underline": "\033[4m"
}
