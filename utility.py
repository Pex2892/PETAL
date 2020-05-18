import globals as gl
import bs4 as bs
import requests
import os
import shutil
import multiprocessing as mlp
import gzip
import hashlib
import configparser
from dateutil.parser import parse
from progressbar import Bar, Counter, ETA, Percentage, ProgressBar
from xml.dom.minidom import parseString
from datetime import datetime


def read_config():
    config = configparser.ConfigParser()
    config.read(gl.filename_config)

    gl.num_cores_input = config['analysis'].getint('n_cpu')
    if gl.num_cores_input == 0 or gl.num_cores_input > mlp.cpu_count():
        # print('Selected all CPUs or an excessive number')
        gl.num_cores_input = mlp.cpu_count()
    print("#CPUs: %d" % gl.num_cores_input)

    gl.pathway_input = config['analysis'].get('pathway')
    print("Pathway: %s" % gl.pathway_input)

    gl.gene_input = config['analysis'].get('gene')
    print("Gene: %s" % gl.gene_input)

    gl.deep_input = config['analysis'].getint('deep')
    print("Deep: %s" % gl.deep_input)


def clear_previous_results():
    pathdir = os.path.join(os.getcwd(), 'export_data')

    if os.path.exists(pathdir):
        shutil.rmtree(pathdir)

    os.makedirs(pathdir)


# da fare
def check_pathway_update_history(url):
    # download and return list of the updated pathway

    pathfile = os.path.join(os.getcwd(), 'database')
    filename = 'pathway_update_history.html.gz'

    # If it does not exist, download the html file with all the pathways recently updated
    download_file(url, pathfile, filename)

    # retrieve datetime to file
    filetime = os.path.getmtime(os.path.join(pathfile, filename))

    # check different time (in seconds) > 24h (seconds), download file
    if (datetime.now() - datetime.fromtimestamp(filetime)).total_seconds() > 86400:
        print('Download updated list, because the saved file has not been updated for more than 24 hours')

        # deleting the oldest file
        os.remove(os.path.join(pathfile, filename))

        # download the html file with all the pathways recently updated
        download_file(url, pathfile, filename)

    # reading the compressed file
    with gzip.open(os.path.join(pathfile, filename), "rb") as f:
        content = f.read().decode('utf-8')

        soup = bs.BeautifulSoup(content, 'html.parser')

        items = soup.findAll('td')

        i = 0
        while i < len(items):
            # check if the field is a date and if from 2020 onwards
            if is_date(items[i].text) and int(items[i].text[0:4]) >= 2020:
                if 'Deleted; ' in items[i + 3].text:

                    # se esiste, rimuovere il pathway eliminato
                    if os.path.exists(os.path.join(pathfile, 'pathways', 'xml', 'hsa' + items[i + 1].text + '.xml.gz')):
                        print('cancello il pathway rimosso da kegg: %s' % items[i + 1].text)
                        os.remove(os.path.join(pathfile, 'pathways', 'xml', 'hsa' + items[i + 1].text + '.xml.gz'))

                    # se esiste, rimuovere il vecchio pathway con cui è stato unito quello cancellato e scarico
                    # quello aggiornato
                    merged_pathway = items[i + 3].text.split('merged into ')[1]
                    if os.path.exists(os.path.join(pathfile, 'pathways', 'xml', 'hsa' + merged_pathway + '.xml.gz')):
                        # retrieve datetime to file
                        filetime = datetime.fromtimestamp(os.path.getmtime(
                            os.path.join(pathfile, 'pathways', 'xml', 'hsa' + items[i + 1].text + '.xml.gz')))

                        # convert time
                        update_time_from_kegg = datetime.strptime(
                            items[i].text + ' 00:00:00.000000', '%Y-%m-%d %H:%M:%S.%f')

                        if filetime < update_time_from_kegg:
                            print('il pathway di kegg è più aggiornato')

                            print(
                                'cancello il pathway vecchio (locale), ma aggiornato su kegg: %s' % 'hsa' + merged_pathway)
                            os.remove(os.path.join(pathfile, 'pathways', 'xml', 'hsa' + merged_pathway + '.xml.gz'))

                            # scaricare il pathway con cui è stato unito
                            print('scarico il pathway aggiornato da kegg: %s' % 'hsa' + merged_pathway)
                            download_file('http://rest.kegg.jp/get/' + 'hsa' + merged_pathway + '/kgml',
                                          os.path.join(pathfile, 'pathways', 'xml'),
                                          'hsa' + merged_pathway + '.xml.gz')
                        else:
                            print('il pathway locale è più aggiornato')

                elif 'Newly added' in items[i + 3].text:
                    print('aggiunto')

                    # verificare se esiste il file
                    if os.path.exists(os.path.join(pathfile, 'pathways', 'xml', 'hsa' + items[i + 1].text + '.xml.gz')):
                        # retrieve datetime to file
                        filetime = datetime.fromtimestamp(os.path.getmtime(
                            os.path.join(pathfile, 'pathways', 'xml', 'hsa' + items[i + 1].text + '.xml.gz')))

                        # convert time
                        update_time_from_kegg = datetime.strptime(
                            items[i].text + ' 00:00:00.000000', '%Y-%m-%d %H:%M:%S.%f')

                        if filetime < update_time_from_kegg:
                            print('il pathway di kegg è più aggiornato')

                            # rimuovere il vecchio pathway
                            print('cancello il pathway vecchio (locale), ma aggiornato su kegg: %s' % 'hsa' + items[
                                i + 1].text)
                            os.remove(os.path.join(pathfile, 'pathways', 'xml', 'hsa' + items[i + 1].text + '.xml.gz'))

                            # scaricare il pathway aggiornato
                            print('scarico il pathway aggiornato da kegg: %s' % 'hsa' + items[i + 1].text)
                            download_file('http://rest.kegg.jp/get/' + 'hsa' + items[i + 1].text + '/kgml',
                                          os.path.join(pathfile, 'pathways', 'xml'),
                                          'hsa' + items[i + 1].text + '.xml.gz')
                        else:
                            print('il pathway locale è più aggiornato')

                i = i + 4
            else:
                break


def download_file(url, pathfile, filename):
    if not os.path.exists(os.path.join(pathfile, filename)):
        try:
            req = requests.get(url)

            with gzip.open(os.path.join(pathfile, filename), "wb") as f:
                f.write(req.content)

        except requests.exceptions.ConnectionError:
            print("ERROR: Connection refused from KEGG")
            exit(1)


def get_info_gene_initial(pathway_hsa, name_gene):
    filename = os.path.join(os.getcwd(), 'database', 'pathways', 'xml', pathway_hsa + '.xml.gz')

    with gzip.open(filename, "rb") as f:
        content = f.read().decode('utf-8')
        mydoc = parseString(content)

        # GENES
        entry = mydoc.getElementsByTagName('entry')

        for elem in entry:
            string_check = name_gene + ','
            if elem.attributes['type'].value == 'gene' and string_check in \
                    elem.getElementsByTagName('graphics')[0].attributes['name'].value:
                return elem.attributes['name'].value, elem.attributes['link'].value

        print('The gene entered not exit into pathway selected!')
        exit(1)


def download_read_html(url):
    filename = url.split('?')[1].replace('+', '_').replace(':', '')
    filename = hashlib.md5(filename.encode("utf-8")).hexdigest() + '.html.gz'

    # download of the html page containing all the genes, passed as parameters in the url
    download_file(url, os.path.join(os.getcwd(), 'database', 'pathways', 'html'), filename)

    # open pathway of the genes in html
    with gzip.open(os.path.join(os.getcwd(), 'database', 'pathways', 'html', filename), "rb") as f:
        content = f.read().decode('utf-8')

        soup = bs.BeautifulSoup(content, 'html.parser')

        list_pathway = list()
        for link in soup.findAll('a'):
            a_tag = str(link.get('href'))
            if "/kegg-bin/show_pathway" in a_tag:
                a_tag = a_tag.split("+")
                a_tag = a_tag[0].split("?")
                list_pathway.append(a_tag[1])

        # generate a unique elements
        list_pathway = list(set(list_pathway))
    return list_pathway


def is_date(string, fuzzy=False):
    """
    Return whether the string can be interpreted as a date.

    :param string: str, string to check for date
    :param fuzzy: bool, ignore unknown tokens in string if True
    """
    try:
        parse(string, fuzzy=fuzzy)
        return True

    except ValueError:
        return False


def search_id_to_hsa(list_genes_this_pathway, hsa_gene):
    return [item for item in list_genes_this_pathway if hsa_gene == item[1]]


def search_gene_to_id(list_genes_this_pathway, id_gene):
    return [item for item in list_genes_this_pathway if id_gene in item]


def concat_multiple_subtype(list_subtype):
    # in concat_multiple_subtype, concateno tutti i subtype della relazione analizzata
    if len(list_subtype) > 0:
        subtype = []
        for item in list_subtype:
            subtype.append(item.attributes['name'].value)
        return '//'.join(subtype)
    return 'None'


def read_kgml(deep, pathway_hsa, name_gene_start, hsa_gene_start, path, occu):
    filename = os.path.join(os.getcwd(), 'database', 'pathways', 'xml', pathway_hsa + '.xml.gz')

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
                # salvo i geni che hanno hsa e sono di tipo gene
                list_genes_this_pathway.append((elem.attributes['id'].value, elem.attributes['name'].value,
                                                elem2.attributes['name'].value.split(',')[0],
                                                elem.attributes['link'].value))

        # cerco gli id all'interno della mappa per quello specifico hsa(gene)
        # perchè in ogni pathway è diverso, anche se è lo stesso gene
        # se li trova restituisce tuple di diversi id che fanno riferimento al gene
        list_ids_gene_input = search_id_to_hsa(list_genes_this_pathway, hsa_gene_start)

        list_rows = list()
        if len(list_ids_gene_input) > 0:

            # scorro tutte le relazioni all'interno della mappa
            for elem in relation:
                # scorro tutta la lista degli id che fanno riferimento allo stesso gene della stessa mappa
                for id_gene in list_ids_gene_input:
                    # verifico se l'entry2(id di partenza) è uguale a uno degli id salvati nella lista
                    if elem.attributes['entry2'].value == id_gene[0]:

                        # estraggo le info dei geni che ha una relazione con id_gene tornato vuoto perchè potrebbe
                        # essere che l'entry1 non esiste o non è di tipo gene (group, compund)
                        list_gene_relation = search_gene_to_id(list_genes_this_pathway, elem.attributes['entry1'].value)

                        # verifico se esiste una relazione effettivamente
                        if len(list_gene_relation) > 0:
                            row = {
                                'deep': deep,
                                'name_father': name_gene_start,
                                'hsa_father': hsa_gene_start,
                                'name_son': list_gene_relation[0][2],
                                'hsa_son': list_gene_relation[0][1],
                                'url_kegg_son': list_gene_relation[0][3],
                                'relation': elem.attributes['type'].value,
                                'type_rel': concat_multiple_subtype(elem.getElementsByTagName('subtype')),
                                'pathway_of_origin': pathway_hsa,
                                'fullpath': path + '/' + list_gene_relation[0][2],
                                'occurrences': occu
                            }
                            list_rows.append(row)
        return list_rows


def analysis_deep_n(deep, gene, gene_hsa, pathway_this_gene, path, occu):
    download_file('http://rest.kegg.jp/get/' + pathway_this_gene + '/kgml',
                  os.path.join(os.getcwd(), 'database', 'pathways', 'xml'), pathway_this_gene + '.xml.gz')

    list_rows = read_kgml(deep, pathway_this_gene, gene, gene_hsa, path, occu)

    return list_rows


def unified(list_rows_returned):
    # Add the results obtained by the threads in parallel in the main data frame
    for row in list_rows_returned:
        if row is not None:
            for cell in row:
                gl.DF_TREE = gl.DF_TREE.append(cell, ignore_index=True)


def get_info_row_duplicated(df_filtered, gene):
    grouped_df = (df_filtered[df_filtered['name_son'] == gene]).groupby('name_father')

    list_to_do_df = list()
    for key, group in grouped_df:
        hsa_end_refactor = ' '.join(group['hsa_son'].tolist())
        hsa_end_refactor = sorted(set(hsa_end_refactor.split(' ')))
        hsa_end_refactor = ' '.join(hsa_end_refactor)

        url_gene_end_refactor = 'https://www.kegg.jp/dbget-bin/www_bget?' + hsa_end_refactor.replace(' ', '+')
        occurences_calculated = group.iloc[0]['occurrences'] * group.shape[0]

        # riga da aggiornare (prendo la prima perchè devo aggiornare tutti i campi)
        # lista di righe da rimuovere meno quella da conservare
        # nuovo hsa_end del gene
        # nuovo url_end del gene
        # unisco tutte le relation in base ai pathway di origine
        # unisco tutte le type_rel in base ai pathway di origine
        # unisco tutti i pathway_origine
        list_to_do_df.append((
                group.index[0],
                list(filter(group.index[0].__ne__, group.index.values.tolist())),
                hsa_end_refactor,
                url_gene_end_refactor,
                '§§'.join(group['relation'].tolist()),
                '§§'.join(group['type_rel'].tolist()),
                '§§'.join(group['pathway_of_origin'].tolist()),
                occurences_calculated
            )
        )
    return list_to_do_df


def clean_update_row_duplicates(list_to_do_df):
    for row in list_to_do_df:
        for cell in row:
            # aggiorno la riga selezionata, con le nuove stringhe nelle 4 colonne
            cols = ['hsa_son', 'url_kegg_son', 'relation', 'type_rel', 'pathway_of_origin', 'occurrences']

            gl.DF_TREE.loc[cell[0], cols] = [cell[2], cell[3], cell[4], cell[5], cell[6], cell[7]]

            # rimuovo le righe selezionate
            if len(cell[1]) > 0:
                gl.DF_TREE = gl.DF_TREE.drop(cell[1])


def set_progress_bar(action, max_elem):
    """
    Progressbar can't guess maxval.

    :param action: str, ....
    :param max_elem: str, ......
    """
    pb = ProgressBar(widgets=[action, ' ', Percentage(), ' (', Counter(), ' of ',
                              max_elem, ') ', Bar('#'), ' ', ETA()], maxval=int(max_elem))
    return pb
