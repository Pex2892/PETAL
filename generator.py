import pandas as pd
import os
import json


COLUMNS_DF = ['hop', 'name_start', 'hsa_start', 'name_end', 'hsa_end', 'url_hsa_end', 'relation',
              'type_rel', 'pathway_origin']

def init_generate_output(hop):
    # carico il csv del livello attuale con i duplicati
    df_1 = pd.read_csv('results/execution/results_level1.csv', sep=';', header=None, names=COLUMNS_DF)
    df_1 = df_1.sort_values(by=['name_end'])

    len_df_1 = df_1.shape[0]

    lista = [[]]
    for i in range(0, len_df_1):
        lista[0].append(df_1.iloc[i, 1] + ';;' + df_1.iloc[i, 3] + '||' + df_1.iloc[i, 5] + '#' + df_1.iloc[i, 6]
                        + '#--' + df_1.iloc[i, 8])

    # print('Livello 1')
    # print(lista)
    # print(len(lista[0]))
    # print('---------------')

    list_returned = []
    for i in range(2, hop+1):
        # print('Livello ' + str(i))
        list_returned = process_level(i, lista[i - 2])
        # print(list_returned)
        # print(len(list_returned))
        lista.insert(i - 1, list_returned)

    # converto le liste innestate, in una lista
    output = ['\n'.join([str(c) for c in lst]) for lst in lista]

    # converto la lista in una stringa
    output = '\n'.join(output)

    # scrivo il file di output
    with open('results/execution/output_text.txt', 'w') as f:
        f.write(output)


# ----------------------


def process_level(level, list_level_prev):
    # carico il csv del livello attuale
    df = pd.read_csv(os.path.join('results', 'execution', 'results_level'+str(level)+'.csv'), sep=';', header=None, names=COLUMNS_DF)

    # ordino il df per nome gene iniziale
    df = df.sort_values(by=['name_start'])

    # numero di righe del livello attuale
    len_df = df.shape[0]

    # numero di geni del livello precedente
    len_list_level_prev = len(list_level_prev)

    result_list_actual = []

    for i in range(0, len_df):
        gene = ';;' + df.iloc[i, 1] + '||'
        for j in range(0, len_list_level_prev):
            if gene in list_level_prev[j]:
                row = list_level_prev[j].split('||')[0]
                row = row + ';;' + df.iloc[i, 3] + '||' + df.iloc[i, 5] + '#' + df.iloc[i, 6] + '#--' + df.iloc[i, 8]
                result_list_actual.append(row)

    return result_list_actual


# ----------------------


def check_child(lista, child, row_act, links_):
    index_level = 0
    for i in range(0, len(lista)):
        if lista[i]["nodeName"] == child:
            return lista[i]["children"]
        index_level = i+1

    type_link = links_[row_act].split("#", 1)
    json_child = {}
    json_child["nodeName"] = child
    json_child["code"] = "lorem"
    json_child["version"] = "lorem"
    json_child["name"] = ""
    json_child["label"] = ""
    json_child["type"] = "type1"
    json_child["link"] = {}
    json_child["link"]['name'] = type_link[0]

    type_link[1] = type_link[1].split('--')[0]

    json_child["link"]['nodeName'] = type_link[1].replace("#", " - ", type_link[1].count('#')-1)[:-1]
    json_child["link"]['direction'] = "ASYN"
    json_child["link"]['from'] = links_[row_act].split("--")[1]
    json_child["children"] = []

    lista.append(json_child)

    lista = lista[index_level]["children"]

    return lista

def output_json():
    rows = ''
    links = []
    with open('results/execution/output_text.txt') as f:
        rows = f.readlines()

    for i in range(0, len(rows)):
        part_line = rows[i].split("||")
        links.append(part_line[1].replace("\n", ""))

        rows[i] = part_line[0].split(";;")
        if rows[i][-1] == "":
            rows[i].pop()

    json_dict = {}

    json_root = {}
    json_root["nodeName"] = rows[0][0]
    json_root["code"] = "lorem"
    json_root["version"] = "lorem"
    json_root["name"] = "null"
    json_root["label"] = "null"
    json_root["type"] = "type4"
    # json_root["link"] = {}
    # json_root["link"]['name'] = "Link node 1 to 2.1"
    # json_root["link"]['nodeName'] = "NODE NAME 2.1"
    # json_root["link"]['direction'] = "ASYN"
    json_root["children"] = []

    p_root = json_root["children"]
    genes = list()
    genes.append(rows[0][0])

    for i in range(0, len(rows)):  # scorro tutte le righe i-sime
        for j in range(1, len(rows[i])):
            genes.append(rows[i][j])
            p_root = check_child(p_root, rows[i][j], i, links)

        p_root = json_root["children"]

    json_dict['tree'] = json_root

    with open("demo/flare.json", "w") as outfile:
        json.dump(json_dict, outfile, indent=2, sort_keys=False)