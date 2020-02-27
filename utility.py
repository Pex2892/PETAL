import argparse
import re
from datetime import datetime
import bs4 as bs
import requests
import os
import shutil
from xml.dom import minidom
from Bcolors import Bcolors as Bc



def command_line():
    parser = argparse.ArgumentParser(description='PETAL (Parallel Pathway Analyzer) \n Scrivere una breve descrizione')

    parser.add_argument("-pathway", default='hsa04010', type=str,
                        help="This is the name of the biological pathway", required=True)

    parser.add_argument("-gene", default='MAPK1', type=str,
                        help="This is the name of the gene", required=True)

    parser.add_argument("-hop", default=2, type=int,
                        help="Represents the maximum search depth", required=True)

    args = parser.parse_args()

    return args.pathway, args.gene, args.hop


def clear_results():
    if os.path.isdir('results'):
        shutil.rmtree('results/execution')
        os.makedirs('results/execution')


def get_info_gene_initial(pathway_hsa, name_gene):
    mydoc = minidom.parse('results/pathways_xml/' + pathway_hsa + '.xml')

    # GENES
    entry = mydoc.getElementsByTagName('entry')

    for elem in entry:
        string_check = name_gene + ','
        if elem.attributes['type'].value == 'gene' and string_check in elem.getElementsByTagName('graphics')[0].attributes['name'].value:
            return elem.attributes['name'].value, elem.attributes['link'].value

    print(Bc.FAIL + Bc.BOLD + "The gene entered not exit into pathway selected! " + Bc.ENDC)
    exit(1)


def download_update_pathway_html(url_page):
    # retrieve datetime to file
    file_time = os.path.getmtime('results/pathway_map_update_history.html')

    # check datetime (in seconds) > 24h (seconds)
    if (datetime.now() - datetime.fromtimestamp(file_time)).total_seconds() > 86400:
        # download file html updated
        try:
            req = requests.get(url_page)

            with open('results/pathway_map_update_history.html', 'wb') as f:
                f.write(req.content)

            soup = bs.BeautifulSoup(req.content, 'html.parser')

            return [item.text for item in soup.findAll('td')[1::4]]
        except requests.exceptions.ConnectionError:
            print("Connection refused from KEGG")
            exit(1)
    else:
        # open file html previously downloaded
        with open('results/pathway_map_update_history.html', 'r') as f:
            soup = bs.BeautifulSoup(f.read(), 'html.parser')

            return [item.text for item in soup.findAll('td')[1::4]]

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
    if not os.path.exists('results/pathways_xml/'+pathway_hsa_input+'.xml'):
        # print('file not exist!')
        pathway_kgml = requests.get('http://rest.kegg.jp/get/' + pathway_hsa_input + '/kgml')
        with open('results/pathways_xml/'+pathway_hsa_input+'.xml', 'wb') as f:
            f.write(pathway_kgml.content)

def download_read_html(url_pathway):
    hsa_pathwway = url_pathway.split('?')[1].replace('+', '_').replace(':', '')

    if not os.path.exists('results/pathways_html/'+hsa_pathwway+'.html'):
        try:
            req = requests.get(url_pathway)

            with open('results/pathways_html/'+hsa_pathwway+'.html', 'wb') as f:
                f.write(req.content)
        except requests.exceptions.ConnectionError:
            print("Connection refused")
            exit(1)

    # open pathway of the genes in html
    with open('results/pathways_html/'+hsa_pathwway+'.html', 'r') as f:
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


def search_id_to_hsa(list_genes_pathway, hsa_gene):
    # print('search for name:', name_gene)
    return [item for item in list_genes_pathway if hsa_gene in item[1]]


def search_gene_to_id(list_genes_pathway, id_gene):
    # print('search for id:', id_gene)
    return [item for item in list_genes_pathway if id_gene in item]


def concat_multiple_subtype(list_subtype):
    if len(list_subtype) > 0:
        subtype = []
        for item in list_subtype:
            subtype.append(item.attributes['name'].value)
        return '//'.join(subtype)
    return 'None'


def read_kgml(level, pathway_hsa, name_gene_start, hsa_gene_start):
    mydoc = minidom.parse('results/pathways_xml/' + pathway_hsa + '.xml')

    # GENES
    entry = mydoc.getElementsByTagName('entry')
    graphics = mydoc.getElementsByTagName('graphics')

    # RELATIONS
    relation = mydoc.getElementsByTagName('relation')
    # subtype = mydoc.getElementsByTagName('subtype')

    list_genes_pathway = []
    for elem, elem2 in zip(entry, graphics):
        if 'hsa' in elem.attributes['name'].value and not 'path:' in elem.attributes['name'].value and elem.attributes[
            'type'].value == 'gene':
            list_genes_pathway.append((elem.attributes['id'].value,
                                       elem.attributes['name'].value,
                                       elem2.attributes['name'].value.split(',')[0],
                                       elem.attributes['link'].value))
    # print('LEN list_genes_pathway:', len(list_genes_pathway))
    # print('list_genes_pathway:', list_genes_pathway)

    list_gene_input = search_id_to_hsa(list_genes_pathway, hsa_gene_start)

    # print('list_gene_input:', list_gene_input)
    # print('LEN list_gene_input:', len(list_gene_input))

    if len(list_gene_input) > 0:
        str_to_csv = ''

        for elem in relation:
            for index_gene in list_gene_input:
                if elem.attributes['entry2'].value == index_gene[0]:
                    # print('ID trovato: ' + elem.attributes['entry1'].value)

                    list_gene_relation = search_gene_to_id(list_genes_pathway, elem.attributes['entry1'].value)

                    if len(list_gene_relation) > 0:
                        # nell riga 71 concateno tutti i subtype della relazione analizzata

                        str_to_csv = str_to_csv + str(level) + ';' + \
                                     name_gene_start + ';' + \
                                     hsa_gene_start + ';' + \
                                     list_gene_relation[0][2] + ';' + \
                                     list_gene_relation[0][1] + ';' + \
                                     elem.attributes['type'].value + ';' + \
                                     concat_multiple_subtype(elem.getElementsByTagName('subtype')) + ';' + \
                                     list_gene_relation[0][3] + ';' + pathway_hsa + "\n"

        tfile = open('results/execution/results_level' + str(level) + '.csv', 'a')
        tfile.write(str_to_csv)
        tfile.close()


def execute_i(url_pathway_kegg, pathway_hsa, level_actual, name_gene_start, hsa_gene_start, list_pathways_gene_actual):
    #global list_child_pathway

    # rimuovo dalla lista trovata il pathway in input
    if list_pathways_gene_actual == None:
        #print('EXECUTE_I:', url_pathway_kegg, pathway_hsa, level_actual, name_gene_start, hsa_gene_start)

        list_pathways_gene_actual = [x for x in download_read_html(url_pathway_kegg) if x != pathway_hsa]
        #print('list_pathways_gene_actual_n:', list_pathways_gene_actual)
        for val in list_pathways_gene_actual:
            download_xml(val)
            read_kgml(level_actual, val, name_gene_start, hsa_gene_start)
    else:
        download_xml(list_pathways_gene_actual)
        read_kgml(level_actual, list_pathways_gene_actual, name_gene_start, hsa_gene_start)
        #print('list_pathways_gene_actual_1', list_pathways_gene_actual)

    #list_child_pathway = []