import globals as gl
import bs4 as bs
import requests
import os
import shutil
import multiprocessing as mlp
import gzip
import configparser
import glob
import re
import json
import numpy as np
from progressbar import Bar, Counter, ETA, Percentage, ProgressBar
from pandas import read_csv
from datetime import datetime
from pywget import wget
from zipfile import ZipFile


def read_config():
    config = configparser.ConfigParser()
    config.read(os.path.join(os.getcwd(), gl.filename_config))

    gl.pathway_input = config['analysis'].get('pathway')
    print(f'Pathway: {gl.pathway_input}')

    gl.gene_input = config['analysis'].get('gene')
    print(f'Gene: {gl.gene_input}')

    gl.deep_input = config['analysis'].getint('deep')
    print(f'Depth: {gl.deep_input}')

    gl.num_cores_input = config['analysis'].getint('n_cpu')
    if gl.num_cores_input == 0 or gl.num_cores_input > mlp.cpu_count():
        # print('Selected all CPUs or an excessive number')
        gl.num_cores_input = mlp.cpu_count()
    print(f'#CPUs: {gl.num_cores_input}')

    gl.mode_input = config['analysis'].getboolean('mode')


def clear_previous_results():
    path = os.path.join(os.getcwd(), 'export_data')

    if os.path.exists(path):
        shutil.rmtree(path)

    os.makedirs(path)

    filename = f'{gl.pathway_input}_{gl.gene_input}_{gl.deep_input}.zip'
    if os.path.isfile(os.path.join(os.getcwd(), filename)):
        os.remove(os.path.join(os.getcwd(), filename))


def check_database():
    path = os.path.join(os.getcwd(), 'database')

    # Check if exist the database directory
    if not os.path.exists(path):
        # Check if the database is empty
        print("---> This is the first run of PETAL, so download the database from Github!")

        # The database zip is downloaded from Github
        print('---> The zip download has started!')
        wget.download(gl.url_download_database, os.getcwd())
        print("---> The zip download is complete!")

        # unzip
        with ZipFile(os.path.join(os.getcwd(), 'only_database.zip'), 'r') as zf:
            zf.extractall(os.getcwd())
        print("---> The database is ready for use")

        os.remove(os.path.join(os.getcwd(), 'only_database.zip'))

    print("----- CHECK UPDATED PATHWAYS -----")

    db_info = json.load(open(os.path.join(os.getcwd(), 'database', 'db_info.json')))

    delta_time = (datetime.now() - datetime.strptime(db_info['updated_at'], '%Y-%m-%d %H:%M:%S.%f')).seconds

    # Check if it has been more than 48 hours (172800 seconds) since the last check on the KEGG database
    if delta_time > 172800:
        print('---> It\'s been more than 48 hours since the last check!')

        # Check for updated pathways
        check_history_pathways(db_info['updated_at'])

        # Update the db_info.json file
        update_info_db(db_info['created_at'])

    # Load the complete list of genes (homo sapiens) into memory
    read_list_homo_sapiens_genes()


def check_history_pathways(last_update):
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
                    print(f'---> The "hsa{items[i + 1].text}" pathway has been added to the KEGG database')

                    # check if the file exists locally
                    if os.path.exists(os.path.join(path, 'pathways', 'kgml', f'hsa{items[i + 1].text}.xml.gz')):
                        os.remove(os.path.join(path, 'pathways', 'kgml', f'hsa{items[i + 1].text}.xml.gz'))
                        print(f'------> The "hsa{items[i + 1].text}.xml.gz" (older) file has been deleted!')

                    flag = download_file(f'http://rest.kegg.jp/get/hsa{items[i + 1].text}/kgml',
                                         f'hsa{items[i + 1].text}.xml.gz', os.path.join(path, 'pathways', 'kgml'))
                    if flag:
                        print(f'------> The "hsa{items[i + 1].text}.xml.gz" file has been downloaded!')

                elif 'Deleted; ' in type_update:
                    merged_into_pathway = items[i + 3].text.split('merged into ')[1]

                    print(f'---> The "hsa{items[i + 1].text}" pathway has been removed from the KEGG and merged '
                          f'into the "hsa{merged_into_pathway}" pathway')

                    # check if the file exists locally
                    if os.path.exists(os.path.join(path, 'pathways', 'kgml', f'hsa{items[i + 1].text}.xml.gz')):
                        os.remove(os.path.join(path, 'pathways', 'kgml', f'hsa{items[i + 1].text}.xml.gz'))
                        print(f'------> The "hsa{items[i + 1].text}.xml.gz" file has been deleted!')

                    # check if the "merged_pathway" file exists locally
                    if os.path.exists(os.path.join(path, 'pathways', 'kgml', f'hsa{merged_into_pathway}.xml.gz')):
                        os.remove(os.path.join(path, 'pathways', 'kgml', f'hsa{merged_into_pathway}.xml.gz'))
                        print(f'------> The "hsa{merged_into_pathway}.xml.gz" file has been deleted!')

                    flag = download_file(f'http://rest.kegg.jp/get/hsa{merged_into_pathway}/kgml',
                                         f'hsa{merged_into_pathway}.xml.gz', os.path.join(path, 'pathways', 'kgml'))
                    if flag:
                        print(f'------> The "hsa{merged_into_pathway}.xml.gz" file has been downloaded!')
            i = i + 4


def read_list_homo_sapiens_genes():
    path = os.path.join(os.getcwd(), 'database', 'genes', 'homo_sapiens_genes.csv')

    # convert array numpy 2d in 1d with flatten
    gl.CSV_GENE_HSA = (read_csv(path, sep="\t").to_numpy()).flatten()

    print("----- LOADED LIST OF HUMAN GENES -----")


def update_info_db(_created_at):
    time_now = str(datetime.now())
    data = {'created_at': _created_at, 'updated_at': time_now}

    with open(os.path.join(os.getcwd(), 'database', 'db_info.json'), 'w') as f:
        json.dump(data, f)

    print('---> The "db_info.json" file has been updated with the latest updates!')


def download_file(url, filename, destination, message='ERROR: Connection refused'):
    try:
        req = requests.get(url)

        if req.text:
            with gzip.open(os.path.join(destination, filename), "wb") as f:
                f.write(req.content)
            return True
        else:
            print(f"The url \"{url}\" returned empty content!")
            return False
    except requests.exceptions.ConnectionError:
        print(message)
        exit(1)


def read_gene_txt(hsa):
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


def API_KEGG_get_name_gene_from_hsa(list_hsa, csv_gene_hsa):
    list_genes_finded = []
    for i_hsa in list_hsa:
        index = np.where(csv_gene_hsa == i_hsa)[0][0]
        list_genes_finded.append(f'{csv_gene_hsa[index + 1].split(";", 1)[0].split(",", 1)[0]}({i_hsa})')

    return ', '.join(list_genes_finded)


def API_KEGG_get_hsa_gene_from_name(gene, csv_gene_hsa):
    r = re.compile(r"^%s[,;]" % gene)
    finded = list(filter(r.match, csv_gene_hsa))[0]
    index = np.where(csv_gene_hsa == finded)[0][0]
    return csv_gene_hsa[index - 1]


def set_progress_bar(action, max_elem):
    pb = ProgressBar(widgets=[action, ' ', Percentage(), ' (', Counter(), ' of ',
                              max_elem, ') ', Bar('#'), ' ', ETA()], maxval=int(max_elem))
    return pb


def export_data_for_deep(deep):
    gl.DF_TREE.to_csv(os.path.join(os.getcwd(), 'export_data', f'df_resulted_deep_{deep}.csv'), sep=';',
                      header=False, index=False)

    gl.DF_TREE.to_csv(os.path.join(os.getcwd(), 'export_data', 'df_resulted.csv'), sep=';', mode='a',
                      header=False, index=False)

    print(f'----- Deep {deep} - CSV SAVED ----- ')


def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


def load_last_csv():
    files = glob.glob(os.path.join(os.getcwd(), 'export_data', 'df_resulted_deep_*.csv'))

    path_last_csv = sorted(files, key=numericalSort)[-1]
    filename_last_csv = os.path.basename(path_last_csv)
    print(f'Found last csv: {filename_last_csv}')

    deep_last_csv = int(re.compile(r'\d+').findall(filename_last_csv)[0])

    if deep_last_csv == gl.deep_input:
        print('The maximum depth selected has already been analyzed.\n'
              'You will find the results in the "exporta_data" folder.')
        exit(1)
    elif deep_last_csv > gl.deep_input:
        print('Warning!\nThe maximum depth analyzed is greater than that selected in input.')
        exit(1)
    elif deep_last_csv < gl.deep_input:
        gl.DF_TREE = read_csv(path_last_csv, sep=";", names=gl.COLS_DF)
        print(f"The analysis will start from the depth of {(deep_last_csv + 1)}")
        return deep_last_csv + 1


def create_zip():
    path = os.path.join(os.getcwd(), 'export_data')
    namezip = f'{gl.pathway_input}_{gl.gene_input}_{gl.deep_input}.zip'

    # create a ZipFile object
    with ZipFile(namezip, 'w') as zipObj:
        # Iterate over all the files in directory
        for folderName, subfolders, filenames in os.walk(path):
            for filename in filenames:
                # create complete filepath of file in directory
                filepath = os.path.join(folderName, filename)
                # Add file to zip
                zipObj.write(filepath, os.path.basename(filepath))
