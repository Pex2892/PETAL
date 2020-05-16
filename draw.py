import globals as gl
import json
import os
from joblib import Parallel, delayed
from utility import set_progress_bar


def draw_json_run():
    """
    Questo metodo in partenza definisce il nodo root, vengono filtrati i geni
    per livello e elaborati in parallelo. Infine il dict globale viene salvato
    nel formato json.

    :param: void.
    """
    # creo nodo root nel dict
    gl.json_dict = {
        'name': gl.gene_input,
        'hsa': gl.gene_input_hsa,
        'url': gl.gene_input_url,
        'info': 'Nothing',
        'occurrences': 0,
        'hop': 0,
        'children': []}

    # ciclo per numero di hop in input
    for i in set_progress_bar('[Drawing]', str(gl.hop_input))(range(1, gl.hop_input+1)):
        # filtro i geni per hop
        df_filter = gl.DF_TREE[gl.DF_TREE['deep'] == i]

        # avvio in parallelo l'aggiunta dei geni nel dict del
        Parallel(n_jobs=gl.num_cores_input, backend='threading')(
            delayed(draw_hop_n)(i, item) for key, item in df_filter.iterrows())

    # salvo il dict in json
    with open(os.path.join(os.getcwd(), 'results', 'data-flare.json'), 'w') as outfile:
        json.dump(gl.json_dict, outfile, indent=4)


def draw_hop_n(hop, item):
    """
    Questo metodo, nel pointer restituito si appende il gene passato al dict globale.

    :param item: pandas.Series, contiene tutte le info del gene passato.
    """
    # prendo il pointer dove salvare i geni
    p = search_key(item['fullpath'], gl.json_dict['children'])
    # inserisco i geni in quel punto trovato
    p.append({
        'name': item['name_son'],
        'hsa': item['hsa_son'],
        'url': item['url_kegg_son'],
        'info': concat_info(item['relation'], item['type_rel'], item['pathway_of_origin']),
        'occurrences': item['occurrences'],
        'hop': hop,
        'children': []
    })


def search_key(path, p):
    """
    Questo metodo, partendo dalla root, cerca livello per livello in base alla variabile "path_arr"
    di trovare la posizione del gene (padre) e posizionare il pointer in quel punto.

    Ritorna il puntatore alla lista "children" del gene padre.

    :param path: str, string to rappresent path this gene. Ex: GeneA/GeneB/GeneC.
    :param p: list, parte dalla lista "children" del gene root.
    """
    path_arr = path.split('/')
    for elem_path in path_arr[1:]:
        for i, item in enumerate(p):
            if item['name'] == elem_path:
                p = p[i]['children']
                break
    return p


def concat_info(rel, type_rel, patwhay):
    """
    Questo metodo, concatena le info (relations, type relations and pathway origin).

    Ritorna una stringa con le info contenate.

    :param rel: list, lista che contiene le relazioni
    :param type_rel: list, lista che contiene le tiplogie di relazioni.
    :param patwhay: list, lista che contiene pathway di origine.
    """
    rel_arr = rel.split('§§')
    type_rel_arr = type_rel.split('§§')
    patwhay_arr = patwhay.split('§§')

    str_info = ' - '.join([c + ' | ' + a + ' | ' + b for a, b, c in zip(rel_arr, type_rel_arr, patwhay_arr)])

    return str_info
