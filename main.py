import shutil
import time
from joblib import Parallel, delayed
import multiprocessing
import subprocess
from xml.dom import minidom
import pandas as pd
import os
import bs4 as bs
import requests
from generator import *

start_time = time.time()

# ----- INITIAL CONFIG ----
pathway_hsa = 'hsa04010'
name_gene = 'MAPK1'
hsa_gene = 'hsa:5594'  # non deve essere statico, ma prelevarlo da questo url:
                       # https://www.genome.jp/dbget-bin/www_bget?pathway+hsa04010

# https://www.genome.jp/dbget-bin/www_bget?hsa:5594+hsa:5595
hop = 2


# ----- INITIAL VARIABLE ----
pd.set_option('display.max_colwidth', -1) # mi serve a non troncare le info salvate all'interno del DF
col = ['hop', 'name_start', 'hsa_start', 'name_end', 'hsa_end', 'url_hsa_end', 'relation', 'type_rel', 'pathway_origin']
df_tree = pd.DataFrame(columns=col)

# ----- INITIAL METHODS ----

def download_xml(pathway_hsa_input):
    if not os.path.exists('./pathways_xml/'+pathway_hsa_input+'.xml'):
        # print('file not exist!')
        pathway_kgml = requests.get('http://rest.kegg.jp/get/' + pathway_hsa_input + '/kgml')
        with open('pathways_xml/'+pathway_hsa_input+'.xml', 'wb') as f:
            f.write(pathway_kgml.content)


def read_html(url_pathway):
    sauce = requests.get(url_pathway)
    soup = bs.BeautifulSoup(sauce.content, 'html.parser')

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


def read_kgml(level, pathway_hsa, name_gene_start, hsa_gene_start):
    list_child_pathway = []
    str_to_csv = ''

    mydoc = minidom.parse('pathways_xml/'+pathway_hsa+'.xml')

    # GENESE
    entry = mydoc.getElementsByTagName('entry')
    graphics = mydoc.getElementsByTagName('graphics')

    # RELATIONS
    relation = mydoc.getElementsByTagName('relation')
    subtype = mydoc.getElementsByTagName('subtype')

    list_genes_pathway = []
    for elem, elem2 in zip(entry, graphics):
        if 'hsa' in elem.attributes['name'].value and not 'path:' in elem.attributes['name'].value and elem.attributes['type'].value == 'gene':
            list_genes_pathway.append((elem.attributes['id'].value,
                                       elem.attributes['name'].value,
                                       elem2.attributes['name'].value.split(',')[0],
                                       elem.attributes['link'].value))
    #print('list_genes_pathway:', list_genes_pathway)

    list_gene_input = search_id_to_hsa(list_genes_pathway, hsa_gene_start)

    # print('list_gene_input:', list_gene_input)

    if len(list_gene_input) > 0:

        #list_relations_pathway = []

        for elem,elem2 in zip(relation, subtype):

            if elem.attributes['entry2'].value == list_gene_input[0][0]:

                list_gene_relation = search_gene_to_id(list_genes_pathway, elem.attributes['entry1'].value)

                if len(list_gene_relation) > 0 and not list_gene_relation[0][2] in list_child_pathway: # togliendo la condizione dopo l'and vedrei i duplicati, quindi fare il calcolo delle occorrenze
                    str_to_csv = str_to_csv + str(level)+';'+name_gene_start+';'+hsa_gene_start+';'+list_gene_relation[0][2]+';'+\
                                 list_gene_relation[0][1]+';'+elem.attributes['type'].value+';'+elem2.attributes['name'].value+';'+list_gene_relation[0][3]+';'+pathway_hsa+"\n"

                    list_child_pathway.append(list_gene_relation[0][2])

        #print(str_to_csv)
        tfile = open('results/results_level' + str(level) + '.csv', 'a')
        tfile.write(str_to_csv)
        tfile.close()

def execute_i(url_pathway_kegg, pathway_hsa, level_actual, name_gene_start, hsa_gene_start, list_pathways_gene_actual):
    global list_child_pathway

    # rimuovo dalla lista trovata il pathway in input
    if list_pathways_gene_actual == None:
        #print('EXECUTE_I:', url_pathway_kegg, pathway_hsa, level_actual, name_gene_start, hsa_gene_start)

        list_pathways_gene_actual = [x for x in read_html(url_pathway_kegg) if x != pathway_hsa]
        #print('list_pathways_gene_actual_n:', list_pathways_gene_actual)
        for val in list_pathways_gene_actual:
            download_xml(val)
            read_kgml(level_actual, val, name_gene_start, hsa_gene_start)
    else:
        download_xml(list_pathways_gene_actual)
        read_kgml(level_actual, list_pathways_gene_actual, name_gene_start, hsa_gene_start)
        #print('list_pathways_gene_actual_1', list_pathways_gene_actual)

    list_child_pathway = []


# ----- INITIAL MAIN ----
if os.path.isdir('results'):
    shutil.rmtree('results')
os.makedirs('results')

num_cores = multiprocessing.cpu_count()
for level_actual in range(1, hop+1): #hop+1 perchè non devo saltare l'ultimo livello
    if level_actual == 1:
        download_xml(pathway_hsa)
        read_kgml(level_actual, pathway_hsa, name_gene, hsa_gene)

        list_pathways_gene_actual_first_level = [x for x in read_html('https://www.genome.jp/dbget-bin/www_bget?hsa:5594+hsa:5595') if x != pathway_hsa]

        inputs = range(0, len(list_pathways_gene_actual_first_level))

        Parallel(n_jobs=num_cores)(delayed(execute_i)('https://www.genome.jp/dbget-bin/www_bget?hsa:5594+hsa:5595', pathway_hsa, level_actual, name_gene, hsa_gene, list_pathways_gene_actual_first_level[i]) for i in inputs)

        #print('EXECUTE_1:', pathway_hsa, level_actual, pathway_hsa, name_gene, hsa_gene)

    else:
        df_genes_level_prev = (df_tree[df_tree['hop'] == level_actual - 1])#.drop_duplicates(subset='name_end')

        inputs = range(0, df_genes_level_prev.shape[0])

        #, verbose=5
        Parallel(n_jobs=num_cores)(delayed(execute_i)(df_genes_level_prev.iloc[i,7], pathway_hsa, level_actual, df_genes_level_prev.iloc[i,3], df_genes_level_prev.iloc[i,4], None) for i in inputs)

    # carico il csv del livello attuale con i duplicati
    df_tree = pd.read_csv('results/results_level'+str(level_actual)+'.csv', sep=';', header=None, names=col)

    # elimino i duplicati dal livello attuale e ri-salvo il file
    df_tree = df_tree.drop_duplicates(subset='name_end')
    df_tree.to_csv('results/results_level'+str(level_actual)+'.csv', sep=';', header=False, index=False)


subprocess.call('sh concatenate.sh', shell=True)

# genero il file di output per poi passarlo e creare il json
init_generate_output(hop)


print('Time of Execution:', time.time()-start_time)