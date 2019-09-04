from xml.dom import minidom
import pandas as pd
import os
import bs4 as bs
import requests

pd.set_option('display.max_colwidth', -1)
df_tree_p = pd.DataFrame(
    columns=['hop', 'name_start', 'hsa_start', 'name_end', 'hsa_end', 'url_hsa_end', 'relation', 'type_rel'])
list_child_pathway = []

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
    global df_tree_p

    mydoc = minidom.parse('pathways_xml/'+pathway_hsa+'.xml')

    # GENESE
    entry = mydoc.getElementsByTagName('entry')
    graphics = mydoc.getElementsByTagName('graphics')

    # RELATIONS
    relation = mydoc.getElementsByTagName('relation')
    subtype = mydoc.getElementsByTagName('subtype')

    # ----- START LIST GENE FOR THIS PATHWAY -----
    # costruisco una lista che contiene di geni con relativo id, hsa, name, link
    list_genes_pathway = []
    for elem, elem2 in zip(entry, graphics):
        if 'hsa' in elem.attributes['name'].value and not 'path:' in elem.attributes['name'].value and elem.attributes['type'].value == 'gene':
            list_genes_pathway.append((elem.attributes['id'].value,
                                       elem.attributes['name'].value,
                                       elem2.attributes['name'].value.split(',')[0],
                                       elem.attributes['link'].value))
    print('list_genes_pathway:', list_genes_pathway)
    # ----- END LIST GENE FOR THIS PATHWAY -----

    # ----- START LIST GENE IN INPUT -----
    # creo una lista che contiene
    list_gene_input = search_id_to_hsa(list_genes_pathway, hsa_gene_start)
    print('list_gene_input:', list_gene_input)
    # ----- END LIST GENE IN INPUT -----

    if len(list_gene_input) > 0:

        list_relations_pathway = []
        for elem,elem2 in zip(relation, subtype):
            if elem.attributes['entry2'].value == list_gene_input[0][0]:
                list_gene_relation = search_gene_to_id(list_genes_pathway, elem.attributes['entry1'].value)
                if name_gene_start == 'BRAF':
                    print('list_gene_relation:', list_gene_relation)
                # print('list_gene_relation:', list_gene_relation)
                # print('list_child_pathway:', list_child_pathway)
                if len(list_gene_relation) > 0 and not list_gene_relation[0][2] in list_child_pathway:
                    # print('HSA CHILD:', list_gene_relation[0][1])
                    # print('CHILD:', list_gene_relation[0][2])
                    list_relations_pathway.append(
                        pd.Series([level, name_gene_start, hsa_gene_start, list_gene_relation[0][2],
                                   list_gene_relation[0][1],
                                   list_gene_relation[0][3],
                                   elem.attributes['type'].value,
                                   elem2.attributes['name'].value
                                   ], index=df_tree_p.columns))
                    list_child_pathway.append(list_gene_relation[0][2])

        # print('list_relations_pathway:', list_relations_pathway)
        if len(list_relations_pathway) > 0:
            df_tree_p = df_tree_p.append(list_relations_pathway, ignore_index=True)

    print(df_tree_p)

######################

def execute(url_pathway_kegg, pathway_hsa, level_actual, name_gene_start, hsa_gene_start):
    list_pathways_gene_actual = read_html(url_pathway_kegg)

    for val in list_pathways_gene_actual:
        if val != pathway_hsa:
            download_xml(val)
            read_kgml(level_actual, val, name_gene_start, hsa_gene_start)

    if len(df_tree_p) > 0:
        if name_gene_start == 'BRAF':
            #print(df_tree_p)
            exit(1)
        return df_tree_p
        #df_tree_p.to_csv('results/'+name_gene_start+'_'+str(level_actual)+'.csv', header=False, index=False)