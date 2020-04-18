import globals as gl
import argparse
import re
from datetime import datetime
import bs4 as bs
import requests
import os
import shutil
from xml.dom import minidom
import json
import time
import multiprocessing as mlp


def read_config():
    with open(gl.filename_config) as f:
        data = json.load(f)

        gl.num_cores_input = int(data['n_CPU']['value'])
        if gl.num_cores_input == 0 or gl.num_cores_input > mlp.cpu_count():
            gl.num_cores_input = mlp.cpu_count()

        gl.log_input = int(data['log']['value'])
        gl.pathway_input = data['pathway']['value']
        gl.gene_input = data['gene']['value']
        gl.hop_input = int(data['hop']['value'])

    print(gl.COLORS['yellow'] + gl.COLORS['bold'] + "--- INITIAL SHELL PARAMETERS ---" + gl.COLORS['end_line'])
    print(gl.COLORS['blue'] + gl.COLORS['bold'] + "#CPUs: %s" % gl.num_cores_input + gl.COLORS['end_line'])
    print(gl.COLORS['blue'] + gl.COLORS['bold'] + "Pathway: %s" % gl.pathway_input + gl.COLORS['end_line'])
    print(gl.COLORS['blue'] + gl.COLORS['bold'] + "Gene: %s" % gl.gene_input + gl.COLORS['end_line'])
    print(gl.COLORS['blue'] + gl.COLORS['bold'] + "Hop: %s" % gl.hop_input + gl.COLORS['end_line'])


def clear_previous_results():
    if os.path.isdir('results'):
        shutil.rmtree('results/execution')
        os.makedirs('results/execution')


def download_update_pathway_html(url_page):
    # download and return list of the updated pathway

    print(gl.COLORS['yellow'] + gl.COLORS['bold'] + "--- CHECK UPDATED PATHWAYS ---" + gl.COLORS['end_line'])

    # retrieve datetime to file
    file_time = os.path.getmtime('results/pathway_map_update_history.html')

    # check datetime (in seconds) > 24h (seconds)
    if (datetime.now() - datetime.fromtimestamp(file_time)).total_seconds() > 86400:
        print(gl.COLORS['green'] + "Download updated list!" + gl.COLORS['end_line'])

        # download file html updated
        try:
            req = requests.get(url_page)

            with open('results/pathway_map_update_history.html', 'wb') as f:
                f.write(req.content)

            soup = bs.BeautifulSoup(req.content, 'html.parser')

            gl.LIST_UPDATED_PATHWAY = [item.text for item in soup.findAll('td')[1::4]]
            # return [item.text for item in soup.findAll('td')[1::4]]
        except requests.exceptions.ConnectionError:
            print("Connection refused from KEGG")
            exit(1)
    else:
        print(gl.COLORS['green'] + "No updates found!" + gl.COLORS['end_line'])

        # open file html previously downloaded
        with open('results/pathway_map_update_history.html', 'r') as f:
            soup = bs.BeautifulSoup(f.read(), 'html.parser')

            gl.LIST_UPDATED_PATHWAY = [item.text for item in soup.findAll('td')[1::4]]

            # return [item.text for item in soup.findAll('td')[1::4]]


# def command_line():
#     parser = argparse.ArgumentParser(description='PETAL (Parallel Pathway Analyzer) \n Scrivere una breve descrizione')
#
#     parser.add_argument("-pathway", default='hsa04010', type=str,
#                         help="This is the name of the biological pathway", required=True)
#
#     parser.add_argument("-gene", default='MAPK1', type=str,
#                         help="This is the name of the gene", required=True)
#
#     parser.add_argument("-hop", default=2, type=int,
#                         help="Represents the maximum search depth", required=True)
#
#     args = parser.parse_args()
#
#     return args.pathway, args.gene, args.hop


def get_info_gene_initial(pathway_hsa, name_gene):
    mydoc = minidom.parse('results/pathways_xml/' + pathway_hsa + '.xml')

    # GENES
    entry = mydoc.getElementsByTagName('entry')

    for elem in entry:
        string_check = name_gene + ','
        if elem.attributes['type'].value == 'gene' and string_check in \
                elem.getElementsByTagName('graphics')[0].attributes['name'].value:
            return elem.attributes['name'].value, elem.attributes['link'].value

    print(gl.COLORS['red'] + gl.COLORS['bold'] + "The gene entered not exit into pathway selected! " +
          gl.COLORS['end_line'])
    exit(1)


# def check_update_pathway_html(check_pathway):
# print(check_pathway + ' - ' + re.findall(r'\d+', check_pathway)[0])
# with open('results/pathway_map_update_history.html', 'r') as f:
#     soup = bs.BeautifulSoup(f.read(), 'html.parser')
#
#     for item in soup.findAll('td')[1::4]:
#         if item.text in check_pathway:
#             return True
# return False


def download_xml(pathway_hsa_input):
    # pathway non esiste
    if not os.path.exists('results/pathways_xml/' + pathway_hsa_input + '.xml'):
        # print('file not exist!')
        pathway_kgml = requests.get('http://rest.kegg.jp/get/' + pathway_hsa_input + '/kgml')
        with open('results/pathways_xml/' + pathway_hsa_input + '.xml', 'wb') as f:
            f.write(pathway_kgml.content)


def download_read_html(url_pathway):
    hsa_pathwway = url_pathway.split('?')[1].replace('+', '_').replace(':', '')

    if not os.path.exists('results/pathways_html/' + hsa_pathwway + '.html'):
        try:
            req = requests.get(url_pathway)

            with open('results/pathways_html/' + hsa_pathwway + '.html', 'wb') as f:
                f.write(req.content)
        except requests.exceptions.ConnectionError:
            print("Connection refused")
            exit(1)

    # open pathway of the genes in html
    with open('results/pathways_html/' + hsa_pathwway + '.html', 'r') as f:
        soup = bs.BeautifulSoup(f.read(), 'html.parser')

        list_pathway = list()
        for link in soup.findAll('a'):
            a_tag = str(link.get('href'))
            if "/kegg-bin/show_pathway" in a_tag:
                a_tag = a_tag.split("+")
                a_tag = a_tag[0].split("?")
                list_pathway.append(a_tag[1])

        list_pathway = list(set(list_pathway))  # generate a unique elements

    return list_pathway


def search_id_to_hsa(list_genes_this_pathway, hsa_gene):
    # print('search for name:', name_gene)
    return [item for item in list_genes_this_pathway if hsa_gene == item[1]]


def search_gene_to_id(list_genes_this_pathway, id_gene):
    # print('search for id:', id_gene)
    return [item for item in list_genes_this_pathway if id_gene in item]


def concat_multiple_subtype(list_subtype):
    if len(list_subtype) > 0:
        subtype = []
        for item in list_subtype:
            subtype.append(item.attributes['name'].value)
        return '//'.join(subtype)
    return 'None'


def read_kgml(i, hop, pathway_hsa, name_gene_start, hsa_gene_start):
    # print('[%s][hop: %s][pathway: %s][gene: %s]\n' % (i, hop, pathway_hsa, name_gene_start))

    mydoc = minidom.parse('results/pathways_xml/' + pathway_hsa + '.xml')

    # GENES
    entry = mydoc.getElementsByTagName('entry')
    graphics = mydoc.getElementsByTagName('graphics')

    # RELATIONS
    relation = mydoc.getElementsByTagName('relation')
    # subtype = mydoc.getElementsByTagName('subtype')

    list_genes_this_pathway = []
    for elem, elem2 in zip(entry, graphics):
        if 'hsa:' in elem.attributes['name'].value and \
                elem.attributes['type'].value == 'gene':
            # not 'path:' in elem.attributes['name'].value and \

            # salvo i geni che hanno hsa e sono di tipo gene
            list_genes_this_pathway.append((elem.attributes['id'].value,
                                            elem.attributes['name'].value,
                                            elem2.attributes['name'].value.split(',')[0],
                                            elem.attributes['link'].value))

    # print('[%s][hop: %s][pathway: %s][gene: %s]list_genes_this_pathway: %s\n' % (i, hop, pathway_hsa, name_gene_start, list_genes_this_pathway))

    # cerco gli id all'interno della mappa per quello specifico hsa(gene)
    # perchè in ogni pathway è diverso, anche se è lo stesso gene
    # se li trova restituisce tuple di diversi id che fanno riferimento al gene
    list_ids_gene_input = search_id_to_hsa(list_genes_this_pathway, hsa_gene_start)

    # print('[%s][hop: %s][pathway: %s][gene: %s]list_ids_gene_input: %s\n' % (i, hop, pathway_hsa, name_gene_start, list_ids_gene_input))

    list_rows = list()
    if len(list_ids_gene_input) > 0:
        # str_to_csv = ''

        # scorro tutte le relazioni all'interno della mappa
        for elem in relation:
            # scorro tutta la lista degli id che fanno riferimento allo stesso gene della stessa mappa
            for id_gene in list_ids_gene_input:
                # verifico se l'entry2(id di partenza) è uguale a uno degli id salvati nella lista
                if elem.attributes['entry2'].value == id_gene[0]:

                    # estraggo le info dei geni che ha una relazione con id_gene
                    # tornato vuoto perchè potrebbe essere che l'entry1 non esiste o non è di tipo gene (group, compund)
                    list_gene_relation = search_gene_to_id(list_genes_this_pathway, elem.attributes['entry1'].value)

                    # print('[%s][hop: %s][pathway: %s][gene: %s]list_gene_relation: %s\n' % (
                    # i, hop, pathway_hsa, name_gene_start, list_gene_relation))

                    # ONLY DEBUG
                    # if len(list_gene_relation) == 0:
                    #    print('mi blocco')

                    # print('FIND! @ ID: ' + elem.attributes['entry1'].value + ' - name: ' + list_gene_relation[0][2])

                    # verifico se esiste una relazione effettivamente
                    if len(list_gene_relation) > 0:
                        # nell riga 71 concateno tutti i subtype della relazione analizzata

                        # str_to_csv = str_to_csv + str(level) + ';' + \
                        #              name_gene_start + ';' + \
                        #              hsa_gene_start + ';' + \
                        #              list_gene_relation[0][2] + ';' + \
                        #              list_gene_relation[0][1] + ';' + \
                        #              elem.attributes['type'].value + ';' + \
                        #              concat_multiple_subtype(elem.getElementsByTagName('subtype')) + ';' + \
                        #              list_gene_relation[0][3] + ';' + pathway_hsa + "\n"
                        row = {
                            'hop': hop,
                            'name_start': name_gene_start,
                            'hsa_start': hsa_gene_start,
                            'name_end': list_gene_relation[0][2],
                            'hsa_end': list_gene_relation[0][1],
                            'url_gene_end': list_gene_relation[0][3],
                            'relation': elem.attributes['type'].value,
                            'type_rel': concat_multiple_subtype(elem.getElementsByTagName('subtype')),
                            'pathway_origin': pathway_hsa
                        }
                        list_rows.append(row)
    return list_rows
    # tfile = open('results/execution/results_level' + str(level) + '.csv', 'a')
    # tfile.write(str_to_csv)
    # tfile.close()


# def execute_i(i, url_pathway_kegg, pathway_hsa, level_actual, name_gene_start, hsa_gene_start,
#               list_pathways_gene_actual):
#     # # rimuovo dalla lista trovata il pathway in input
#     # if list_pathways_gene_actual == None:
#     #     # print('EXECUTE_I:', url_pathway_kegg, pathway_hsa, level_actual, name_gene_start, hsa_gene_start)
#     #
#     #     list_pathways_gene_actual = [x for x in download_read_html(url_pathway_kegg) if x != pathway_hsa]
#     #     # print('list_pathways_gene_actual_n:', list_pathways_gene_actual)
#     #     for val in list_pathways_gene_actual:
#     #         download_xml(val)
#     #         read_kgml(i, level_actual, val, name_gene_start, hsa_gene_start)
#     # else:
#     print('----->>>list_pathways_gene_actual: ', list_pathways_gene_actual)
#     download_xml(list_pathways_gene_actual)
#     list_rows = read_kgml(i, level_actual, list_pathways_gene_actual, name_gene_start, hsa_gene_start)
#     exit(1)
#     return list_rows
#     # print('list_pathways_gene_actual_1', list_pathways_gene_actual)

def analysis_hop_n(iteration, hop, gene, gene_hsa, pathway_this_gene):
    # print('[i: %s][hop: %s][gene: %s][hsa: %s][pathway: %s]\n' % (iteration, hop, gene, gene_hsa, pathway_this_gene))

    download_xml(pathway_this_gene)
    list_rows = read_kgml(iteration, hop, pathway_this_gene, gene, gene_hsa)

    return list_rows


def unified(list_rows_returned):
    for row in list_rows_returned:
        if row is not None:
            for cell in row:
                gl.DF_TREE = gl.DF_TREE.append(cell, ignore_index=True)


def get_info_row_duplicated(df_filtered, gene):
    subdf = df_filtered[df_filtered['name_end'] == gene]

    # riga duplicate ma che ha più info rispetta alle altre
    index_gene_more_info = subdf[subdf['hsa_end'] == subdf['hsa_end'].max()].index[0]

    list_to_do_df = [
        index_gene_more_info,  # riga da aggiornare e conservare
        list(filter(index_gene_more_info.__ne__, subdf.index.values.tolist())),  # lista di righe da rimuovere meno
        # quella da conservare
        '§§'.join(subdf['relation'].tolist()),  # unisco tutte le relation in base ai pathway di origine
        '§§'.join(subdf['type_rel'].tolist()),  # unisco tutte le type_rel in base ai pathway di origine
        '§§'.join(subdf['pathway_origin'].tolist())  # unisco tutti i pathway_origine
    ]

    return list_to_do_df

def clean_update_row_duplicates(list_to_do_df):
    for row in list_to_do_df:
        # print(row, '\n')
        # print(gl.DF_TREE.loc[row[0], ['relation', 'type_rel', 'pathway_origin']])

        # aggiorno la riga selezionata, con le nuove stringhe nelle 3 colonne
        gl.DF_TREE.loc[row[0], ['relation', 'type_rel', 'pathway_origin']] = [row[2], row[3], row[4]]

        # rimuovo le righe selezionate
        gl.DF_TREE = gl.DF_TREE.drop(row[1])
        # print(gl.DF_TREE.loc[row[0], ['relation', 'type_rel', 'pathway_origin']])
