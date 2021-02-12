import os
import numpy as np
import gzip
import re
import json
import bs4 as bs
from datetime import datetime
from pywget import wget
from zipfile import ZipFile
from xml.dom.minidom import parseString
from utility import download_file


class KEGG:
    def __init__(self, url: str, delta_time: int):
        self.url = url
        self.delta_time = delta_time  # in seconds –  172800 secs = 48h

        self.check_db()
        self.hs_genes_list = self.read_homo_sapiens_genes()
        self.len_hs_genes_list = len(self.hs_genes_list)

    def check_db(self):
        path = os.path.join(os.getcwd(), 'database')

        # Check if exist the database directory
        if not os.path.exists(path):
            # Check if the database is empty
            print(">>>>>> This is the first run of PETAL, so download the database from Github!")

            # The database zip is downloaded from Github
            print('>>>>>> The zip download has started!')
            wget.download(self.url, os.getcwd())
            print(">>>>>> The zip download is complete!")

            # unzip
            with ZipFile(os.path.join(os.getcwd(), 'only_database.zip'), 'r') as zf:
                zf.extractall(os.getcwd())
            print(">>>>>> The database is ready for use")

            os.remove(os.path.join(os.getcwd(), 'only_database.zip'))

        print(">>> CHECK UPDATED PATHWAYS <<<")

        db_info = json.load(open(os.path.join(os.getcwd(), 'database', 'db_info.json')))

        delta_time = (datetime.now() - datetime.strptime(db_info['updated_at'], '%Y-%m-%d %H:%M:%S.%f')).total_seconds()

        if delta_time > self.delta_time:
            print(f'>>>>>> It\'s been more than {int(self.delta_time/3600)} hours since the last check!')

            # Check for updated pathways
            self.check_history_pathways(db_info['updated_at'])

            # Update the db_info.json file
            self.update_info_db(db_info['created_at'])

    def check_history_pathways(self, last_update):
        url = 'https://www.genome.jp/kegg/docs/upd_map.html'
        path = os.path.join(os.getcwd(), 'database')
        filename = 'history_pathways.html.gz'

        # download the html file with all the pathways recently updated
        download_file(url, filename, path)

        with gzip.open(os.path.join(path, filename), "rb") as f:
            content = f.read().decode('utf-8')

            soup = bs.BeautifulSoup(content, 'html.parser')
            items = soup.findAll('td')

            i = 0
            while i < len(items):
                date_json = datetime.fromisoformat(last_update).strftime("%Y-%m-%d")
                date_pathways = datetime.strptime(items[i].text, "%Y-%m-%d")

                # Exit the loop, if there are updates older than 2020
                if date_pathways.year < 2020:
                    i = len(items)

                if date_json < date_pathways.strftime("%Y-%m-%d"):
                    type_update = items[i + 3].text

                    if 'Newly added' in type_update:
                        print(f'>>>>>> The "hsa{items[i + 1].text}" pathway has been added to the KEGG database')

                        # check if the file exists locally
                        if os.path.exists(os.path.join(path, 'pathways', 'kgml', f'hsa{items[i + 1].text}.xml.gz')):
                            os.remove(os.path.join(path, 'pathways', 'kgml', f'hsa{items[i + 1].text}.xml.gz'))
                            print(f'>>>>>>>>> The "hsa{items[i + 1].text}.xml.gz" (older) file has been deleted!')

                        flag = download_file(f'http://rest.kegg.jp/get/hsa{items[i + 1].text}/kgml',
                                             f'hsa{items[i + 1].text}.xml.gz', os.path.join(path, 'pathways', 'kgml'))
                        if flag:
                            print(f'>>>>>>>>> The "hsa{items[i + 1].text}.xml.gz" file has been downloaded!')

                    elif 'Deleted; ' in type_update:
                        merged_into_pathway = items[i + 3].text.split('merged into ')[1]

                        print(f'>>>>>> The "hsa{items[i + 1].text}" pathway has been removed from the KEGG and merged '
                              f'into the "hsa{merged_into_pathway}" pathway')

                        # check if the file exists locally
                        if os.path.exists(os.path.join(path, 'pathways', 'kgml', f'hsa{items[i + 1].text}.xml.gz')):
                            os.remove(os.path.join(path, 'pathways', 'kgml', f'hsa{items[i + 1].text}.xml.gz'))
                            print(f'>>>>>>>>> The "hsa{items[i + 1].text}.xml.gz" file has been deleted!')

                        # check if the "merged_pathway" file exists locally
                        if os.path.exists(os.path.join(path, 'pathways', 'kgml', f'hsa{merged_into_pathway}.xml.gz')):
                            os.remove(os.path.join(path, 'pathways', 'kgml', f'hsa{merged_into_pathway}.xml.gz'))
                            print(f'>>>>>>>>> The "hsa{merged_into_pathway}.xml.gz" file has been deleted!')

                        flag = download_file(f'http://rest.kegg.jp/get/hsa{merged_into_pathway}/kgml',
                                             f'hsa{merged_into_pathway}.xml.gz', os.path.join(path, 'pathways', 'kgml'))
                        if flag:
                            print(f'>>>>>>>>> The "hsa{merged_into_pathway}.xml.gz" file has been downloaded!')
                i = i + 4
        os.remove(os.path.join(path, filename))

    def update_info_db(self, _created_at):
        time_now = str(datetime.now())
        data = {'created_at': _created_at, 'updated_at': time_now}

        with open(os.path.join(os.getcwd(), 'database', 'db_info.json'), 'w') as f:
            json.dump(data, f)

        print('>>>>>> The "db_info.json" file has been updated with the latest updates!')

    def read_homo_sapiens_genes(self):
        path = os.path.join(os.getcwd(), 'database', 'genes', 'homo_sapiens_genes.csv')

        # convert array numpy 2d in 1d with flatten
        return np.genfromtxt(path, dtype=np.str, delimiter='\t').flatten()

    def get_gene_from_hsa(self, hsa: str):
        for i in range(0, self.len_hs_genes_list, 2):
            if hsa == self.hs_genes_list[i]:
                return [self.hs_genes_list[i], self.get_alias(self.hs_genes_list[i + 1])]

    def get_gene_from_name(self, name: str, flag: bool = False):
        for i in range(0, self.len_hs_genes_list, 2):
            if name in self.hs_genes_list[i + 1]:
                t = self.get_alias(self.hs_genes_list[i + 1])
                for alias in t:
                    if name == alias:
                        return [self.hs_genes_list[i], t,
                                f'https://www.kegg.jp/dbget-bin/www_bget?{self.hs_genes_list[i]}']
        print(f'>>>>>> The gene "{name}" does not exist. Try to verify it manually.')
        if flag:
            exit(0)


    def get_alias(self, alias_list):
        # NOTE: THERE ARE GENES IN THE LIST WITHOUT NAME
        if ';' in alias_list:
            alias = alias_list.split(";", 1)[0]
            if ',' in alias:
                alias = alias.split(", ")
            return alias
        else:
            print('No aliases were found.')
            exit()

    def check_alias(self, name: str, alias: list):
        if name != alias[1][0]:
            print(f'>>>>>> The gene name entered does not match the name provided by KEGG, '
                  f'i.e., "{alias[1][0]}" and not "{name}"')
            return alias[1][0]
        return name

    def check_exist_gene_in_pathway(self, pathway: str, gene: str):
        with gzip.open(os.path.join(os.getcwd(), 'database', 'pathways', 'kgml', f'{pathway}.xml.gz'), "rb") as f:
            content = f.read().decode('utf-8')

        if gene not in content:
            print(f'>>>>>> The "{gene}" gene was not found within the "{pathway}" pathway.')
            exit()

    def read_gene_txt(self, hsa: str):
        # open pathway of the genes in txt
        with gzip.open(os.path.join(os.getcwd(), 'database', 'genes', 'txt', f'{hsa}.gz.txt'), "rb") as f:
            content = f.read().decode('utf-8')

            # Identify the pathways of the gene passed into the method through the regex
            res = re.findall(r'(hsa\d+)? {2}', content)

            # Empty items are removed
            res = list(filter(None, res))

            # Remove duplicates, if they exist
            res = list(set(res))
        return res

    def read_kgml(self, l: list):
        filename = os.path.join(os.getcwd(), 'database', 'pathways', 'kgml', f'{l[1]}.xml.gz')

        with gzip.open(filename, "rb") as f:
            content = f.read().decode('utf-8')
            mydoc = parseString(content)

            # GENES
            entry = mydoc.getElementsByTagName('entry')
            graphics = mydoc.getElementsByTagName('graphics')

            # RELATIONS
            relation = mydoc.getElementsByTagName('relation')

            list_genes_this_pathway = []
            for elem, elem2 in zip(entry, graphics):
                if 'hsa:' in elem.attributes['name'].value and elem.attributes['type'].value == 'gene':
                    # Genes are saved are of the "gene" type
                    list_genes_this_pathway.append((elem.attributes['id'].value, elem.attributes['name'].value,
                                                    elem2.attributes['name'].value.split(',')[0],
                                                    elem.attributes['link'].value))

            # All the ids contained in the pathway are recovered, because a gene can have
            # different ids in different pathways or in the same pathway
            list_ids_gene_input = self.search_id_to_hsa(list_genes_this_pathway, l[3])

            list_rows = []
            if len(list_ids_gene_input) > 0:

                # You scroll through the relationships found within the pathway
                for elem in relation:
                    # Scroll the list of ids found
                    for id_gene in list_ids_gene_input:
                        # Check that "entry2" has at least one match in the id list
                        if elem.attributes['entry1'].value == id_gene[0]:

                            # We try to find a correspondence between gene and your id. If it returns zero,
                            # it indicates that that gene does not exist or is group or compund type and not gene type
                            list_gene_relation = self.search_gene_to_id(list_genes_this_pathway, elem.attributes['entry2'].value)

                            # Checks if at least one connection has been found
                            if len(list_gene_relation) > 0:

                                # It could happen that one final gene, we have many hsa. will be managed individually
                                split_hsa = list_gene_relation[0][1].split(" ")

                                row = {
                                    'depth': l[0],
                                    's_gene': l[2],
                                    's_gene_hsa': l[3],
                                    'e_gene': list_gene_relation[0][2],
                                    'e_gene_hsa': split_hsa[0],
                                    'e_gene_url': f'https://www.kegg.jp/dbget-bin/www_bget?{split_hsa[0]}',
                                    'isoforms': ','.join(split_hsa[1:]),
                                    'relation': elem.attributes['type'].value,
                                    'type_rel': self.concat_multiple_subtype(elem.getElementsByTagName('subtype')),
                                    'pathway_origin': l[1],
                                    'fullpath': f'{l[4]}/{list_gene_relation[0][2]}',
                                    'occurrences': l[5]
                                }
                                list_rows.append(row)
            return list_rows

    def search_id_to_hsa(self, list_genes_this_pathway: list, hsa_gene: str):
        return [item for item in list_genes_this_pathway if hsa_gene in f'{item[1]} ']

    def search_gene_to_id(self, list_genes_this_pathway: list, id_gene: str):
        return [item for item in list_genes_this_pathway if id_gene in item]

    def concat_multiple_subtype(self, list_subtype: list):
        # All the "subtypes" of the selected relationship are concatenated
        if len(list_subtype) > 0:
            subtype = []
            for item in list_subtype:
                subtype.append(item.attributes['name'].value)
            return '//'.join(subtype)
        return 'None'
