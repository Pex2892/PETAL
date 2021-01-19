import globals as gl
import json
import os
import pandas as pd
from joblib import Parallel, delayed
from utility import get_gene_info_from_name, set_progress_bar


def draw_json_run(gene_input, df_resulted):
    """
    This method defines the root node, after which the genes are
    filtered by level and processed in parallel.
    Finally the output is saved in the json format.

    :param: void.
    """

    df = pd.read_csv(df_resulted, sep=";", names=gl.COLS_DF)

    iter_deep = df['deep'].max()
    gene_info = get_gene_info_from_name(gene_input, gl.CSV_GENE_HSA)

    # The root node is created in the dictionary
    gl.json_dict = {
        'name': gene_input,
        'hsa': gene_info[0],
        'url': gene_info[2],
        'info': 'Nothing',
        'occurrences': 0,
        'deep': 0,
        'children': []}

    for i in set_progress_bar('[Drawing]', str(iter_deep))(range(1, iter_deep+1)):
        # Genes divided by depth are filtered
        df_filter = df[df['deep'] == i]

        # The addition of genes in the dictionary by depth starts in parallel
        Parallel(n_jobs=gl.num_cores_input, backend='threading')(
            delayed(draw_deep_n)(i, item) for key, item in df_filter.iterrows())

    # export the dictionary to json
    with open(os.path.join(os.getcwd(), 'demo', 'data-flare.json'), 'w') as outfile:
        json.dump(gl.json_dict, outfile, indent=4)


def draw_deep_n(deep, item):
    # We look for the pointer, where to save the genes
    p = search_key(item['fullpath'], gl.json_dict['children'])
    # The genes found in the selected pointer are inserted
    p.append({
        'name': item['name_son'],
        'hsa': item['hsa_son'],
        'url': item['url_kegg_son'],
        'info': concat_info(item['relation'], item['type_rel'], item['pathway_of_origin']),
        'occurrences': item['occurrences'],
        'deep': deep,
        'children': []
    })


def search_key(path, p):
    # This method, starting from the root, based on the "path_arr" variable, tries
    # to find the position of the gene (father) and position the pointer at that point.
    # Return the pointer to the "children" list of the parent gene.
    path_arr = path.split('/')
    for elem_path in path_arr[1:]:
        for i, item in enumerate(p):
            if item['name'] == elem_path:
                p = p[i]['children']
                break
    return p


def concat_info(rel, type_rel, patwhay):
    rel_arr = rel.split('§§')
    type_rel_arr = type_rel.split('§§')
    patwhay_arr = patwhay.split('§§')

    str_info = ' - '.join([f'{c} | {a} | {b}' for a, b, c in zip(rel_arr, type_rel_arr, patwhay_arr)])

    return str_info
