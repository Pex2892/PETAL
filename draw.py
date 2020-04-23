import globals as gl
from joblib import Parallel, delayed
import json
import os


def run():
    gl.json_dict = {'name': gl.gene_input, 'children': []}

    for i in range(1, gl.hop_input+1):
        df_filter = gl.DF_TREE[gl.DF_TREE['hop'] == i]

        Parallel(n_jobs=gl.num_cores_input, backend='threading')(
            delayed(draw_hop_n)(item) for key, item in df_filter.iterrows())

    with open(os.path.join(os.getcwd(), 'demo', 'data-flare.json'), 'w') as outfile:
        json.dump(gl.json_dict, outfile, indent=4)


def draw_hop_n(item):
    p = search_key(item['path'], gl.json_dict['children'])
    p.append({'name': item['name_end'], 'children': []})


def search_key(path, p):
    path_arr = path.split('/')
    for elem_path in path_arr[1:]:
        for i, item in enumerate(p):
            if item['name'] == elem_path:
                p = p[i]['children']
                break
    return p
