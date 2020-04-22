import globals as gl
import json
import os

def run():
    generate_tree_json()

def generate_tree_json():
    gl.json_dict = {'name': gl.gene_input, 'children': []}

    df_filter = gl.DF_TREE[gl.DF_TREE['hop'] == 1]

    list_t = [
        [],  # salvo l'indice per poter poi recuperare la lista children di questo gene
        []   # salvo i nomi dei geni
    ]
    for key, value in df_filter.iterrows():
        list_t.append((key, value['name_end']))
        gl.json_dict['children'].append({'name': value['name_end'], 'children': []})

        list_t[0].append(key)
        list_t[1].append(value['name_end'])

    hop2(gl.json_dict['children'], list_t)
    #print(list_t)
    #print(list_t[0])
    #print(list_t[1].index('MAP2K2'))



def hop2(p, list_t):
    df_filter = gl.DF_TREE[gl.DF_TREE['hop'] == 2]

    list_t2 = [
        [],  # salvo l'indice per poter poi recuperare la lista children di questo gene
        []   # salvo i nomi dei geni
    ]
    for key, value in df_filter.iterrows():
        # print(key, value)

        search_index_gene = list_t[1].index(value['name_start'])
        get_index_json_dict = list_t[0][search_index_gene]
        print('hop2 - search_index_gene: ', search_index_gene)
        print('hop2 - get_index_json_dict: ', list_t[0][search_index_gene])

        p[get_index_json_dict]['children'].append({'name': value['name_end'], 'children': []})

    with open(os.path.join(os.getcwd(), 'demo', 'data-flare.json'), 'w') as outfile:
        json.dump(gl.json_dict, outfile, indent=4)