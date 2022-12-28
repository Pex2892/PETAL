from os.path import join as pathjoin, exists as pathexists
from os import getcwd, remove
import os
import pandas as pd
import requests
from zipfile import ZipFile
from pandas import read_csv
from lxml import etree
from io import StringIO

def download_to_file(url, filename):
    r = requests.get(url, allow_redirects=True)
    with open(filename, 'wb') as f:
        f.write(r.content)

class Kegg:
    def __init__(self):
        self.url_database = 'https://github.com/Pex2892/PETAL/releases/download/v1.2/only_database.zip'
        self.df_genes = None
        self.df_pathways = None
        self.df_literature = None

    def run(self):
        self._is_downloaded()
        self._load_genes()
        self._load_pathways()
        self._load_genes_literature()
        print('\t\u274C\u0020\u0020!!!! AGGIUNGERE OTTIMIZZAZIONI DEI TIPI https://medium.com/bigdatarepublic/advanced-pandas-optimize-speed-and-memory-a654b53be6c2')

    def _is_downloaded(self):
        p_database = pathjoin(getcwd(), 'database')

        # Check if exist the database directory
        if not pathexists(p_database):
            # Check if the database is empty
            print("\t\u279C This is the first run of PETAL, so download the database from Github!")

            # The database zip is downloaded from Github
            print('\t\u279C The zip download has started!')
            download_to_file(url=self.url_database, 
                filename=os.path.join(os.getcwd(), 'only_database.zip'))           
            print("\t\u279C The zip download is complete!")

            # unzip
            with ZipFile(pathjoin(getcwd(), 'only_database.zip'), 'r') as zf:
                zf.extractall(getcwd())

            remove(pathjoin(getcwd(), 'only_database.zip'))

        print("\t\u279C The database is ready for use")

    def _load_genes(self):
        self.df_genes = read_csv(pathjoin(getcwd(), 'database', 'genes.csv'), sep='\t')
        self.df_genes.set_index('id', inplace=True)
        print("\t\u279C Genes have been loaded")

    def _load_pathways(self):
        self.df_pathways = read_csv(pathjoin(getcwd(), 'database', 'pathways.csv'), sep='\t')
        self.df_pathways.set_index('ID', inplace=True)
        print("\t\u279C Pathways have been loaded")

    def _load_genes_literature(self):
        print('\t\u274C\u0020\u0020!!!! MANCA IL CARICAMENTO DELLA LETTERATURA')

    def read_kgml(self, pathway_id, gene_id, fullpath):
        # print(pathway_id, gene_id, fullpath)
        content = self.df_pathways.at[pathway_id, 'kgml'][2:-1]
        root = etree.parse(StringIO(content))

        xml_entries = root.findall('entry')
        df_entries = self._filter_entries(xml_entries)
        # Search multiple-id from hsa
        indices_entries = df_entries[df_entries['genes_id'].str.contains(f'{gene_id},')].index

        xml_relations = root.findall('relation')
        df_relations = self._filter_relations(xml_relations, indices_entries)

        results = self._build_results(df_entries, df_relations, pathway_id, fullpath)
        return results

    def _filter_entries(self, xml_entries):
        rows = []
        for i in range(0, len(xml_entries)):
            if xml_entries[i].attrib['type'] == 'gene' and 'hsa:' in xml_entries[i].attrib['name']:
                # print(xml_entries[i].attrib['name'], xml_entries[i].attrib['type'])
                rows.append({
                    'xml_id': xml_entries[i].attrib['id'],
                    'genes_id': xml_entries[i].attrib['name'],
                    'url': xml_entries[i].attrib['link'],
                })
        df_entries = pd.DataFrame(rows)
        df_entries.set_index('xml_id', inplace=True)

        # I fix the strings of the genes_id column, so I can search them uniquely
        df_entries['genes_id'] = df_entries['genes_id'].str.replace(' ', ',')
        df_entries['genes_id'] = df_entries['genes_id'] + ','

        return df_entries

    def _filter_relations(self, xml_relations, indices_entries):
        rows = []
        for i in range(0, len(xml_relations)):
            for j in range(0, len(indices_entries)):
                if xml_relations[i].attrib['entry1'] == indices_entries[j]:
                    subtype = [item.attrib['name'] for item in xml_relations[i].findall('subtype')]
                    # print(xml_relations[i].attrib['entry1'], xml_relations[i].attrib['entry2'], xml_relations[i].findall('subtype'))

                    rows.append({
                        'entry1': indices_entries[j],
                        'entry2': xml_relations[i].attrib['entry2'],
                        'relation': xml_relations[i].attrib['type'],
                        'subtype_relation': ';'.join(subtype)
                    })
        df_relations = pd.DataFrame(rows)

        return df_relations

    def _build_results(self, df_entries, df_relations, pathway_id, fullpath):
        rows = []
        for i in range(0, df_relations.shape[0]):
            if df_relations.at[i, 'entry2'] in df_entries.index:
                # print(df_relations.at[i, 'entry1'])
                entry1 = df_entries.at[df_relations.at[i, 'entry1'], 'genes_id'][:-1].split(',')
                # print(df_relations.at[i, 'entry2'])
                entry2 = df_entries.at[df_relations.at[i, 'entry2'], 'genes_id'][:-1].split(',')

                ending_gene_name = self.df_genes.at[entry2[0], 'name'].split(',', 1)[0]

                rows.append({
                    'depth': None,
                    'starting_gene_id': entry1[0],
                    'starting_gene_name': self.df_genes.at[entry1[0], 'name'].split(',', 1)[0],
                    'starting_isoforms_id': ';'.join(entry1[1:]),
                    'starting_isoforms_name': ';'.join(self.df_genes.loc[entry1[1:], 'name'].apply(lambda x: x.split(', ', 1)[0])),
                    'ending_gene_id': entry2[0],
                    'ending_gene_name': ending_gene_name,
                    'ending_isoforms_id': ';'.join(entry2[1:]),
                    'ending_isoforms_name': ';'.join(self.df_genes.loc[entry2[1:], 'name'].apply(lambda x: x.split(', ', 1)[0])),
                    'relation': df_relations.at[i, 'relation'],
                    'subtype': df_relations.at[i, 'subtype_relation'],
                    'reference_pathway': pathway_id,
                    'fullpath': f"{fullpath}/{ending_gene_name}",
                    'occurrences': 1,
                })
        # print(rows)
        return rows


