import globals as gl
from datetime import datetime
import bs4 as bs
import requests
import os
import shutil
from xml.dom import minidom
import json
import time
import multiprocessing as mlp
import gzip
import hashlib
from dateutil.parser import parse
import logging
import logging.config


def read_config():
    with open(gl.filename_config) as f:
        data = json.load(f)

        set_logger(not bool(data['log']['value']))

        gl.logger.info('Reading of the initial parameters')

        gl.num_cores_input = int(data['n_CPU']['value'])
        if gl.num_cores_input == 0 or gl.num_cores_input > mlp.cpu_count():
            gl.logger.info('Selected all CPUs or an excessive number')
            gl.num_cores_input = mlp.cpu_count()
        gl.logger.debug('Number of CPUs: %d' % gl.num_cores_input)

        gl.pathway_input = data['pathway']['value']
        gl.logger.debug('Pathway selected: "%s"' % gl.pathway_input)

        gl.gene_input = data['gene']['value']
        gl.logger.debug('Gene selected: "%s"' + gl.gene_input + '"')

        gl.hop_input = int(data['hop']['value'])
        gl.logger.debug('Maximum depth selected: %d' % gl.hop_input)

    print(gl.COLORS['yellow'] + gl.COLORS['bold'] + "--- INITIAL SHELL PARAMETERS ---" + gl.COLORS['end_line'])
    print(gl.COLORS['blue'] + gl.COLORS['bold'] + "#CPUs: %s" % gl.num_cores_input + gl.COLORS['end_line'])
    print(gl.COLORS['blue'] + gl.COLORS['bold'] + "Pathway: %s" % gl.pathway_input + gl.COLORS['end_line'])
    print(gl.COLORS['blue'] + gl.COLORS['bold'] + "Gene: %s" % gl.gene_input + gl.COLORS['end_line'])
    print(gl.COLORS['blue'] + gl.COLORS['bold'] + "Hop: %s" % gl.hop_input + gl.COLORS['end_line'])


def set_logger(flag):
    gl.logger = logging.getLogger(__name__)
    gl.logger.setLevel(logging.DEBUG)

    # enabled or disabled
    gl.logger.disabled = flag

    fh = logging.FileHandler('PETAL-debug.log', mode='w')
    fh.setLevel(logging.DEBUG)
    gl.logger.addHandler(fh)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    gl.logger.addHandler(fh)


def clear_previous_results():
    pathdir = os.path.join(os.getcwd(), 'results')

    # Remove the previous results folder
    shutil.rmtree(pathdir)
    gl.logger.warning('Removed the previous results folder')

    # Create the results folder
    os.makedirs(pathdir)
    gl.logger.warning('Created the results folder')


# da fare
def check_pathway_update_history(url):
    # download and return list of the updated pathway

    print(gl.COLORS['yellow'] + gl.COLORS['bold'] + "--- CHECK UPDATED PATHWAYS ---" + gl.COLORS['end_line'])
    gl.logger.info('Check for pathway updates')

    pathfile = os.path.join(os.getcwd(), 'database')
    filename = 'pathway_update_history.html.gz'

    # If it does not exist, download the html file with all the pathways recently updated
    download_file(url, pathfile, filename)

    # retrieve datetime to file
    filetime = os.path.getmtime(os.path.join(pathfile, filename))

    # check different time (in seconds) > 24h (seconds), download file
    if (datetime.now() - datetime.fromtimestamp(filetime)).total_seconds() > 86400:
        print(gl.COLORS['green'] + "Download updated list!" + gl.COLORS['end_line'])
        gl.logger.info('The saved file has not been updated for more than 24 hours')

        # deleting the oldest file
        os.remove(os.path.join(pathfile, filename))

        # download the html file with all the pathways recently updated
        download_file(url, pathfile, filename)

    # reading the compressed file
    with gzip.open(os.path.join(pathfile, filename), "rb") as f:
        content = f.read().decode('utf-8')

        soup = bs.BeautifulSoup(content, 'html.parser')

        items = soup.findAll('td')

        gl.logger.debug('Check the pathways updated from 2020 onwards')

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


def get_info_gene_initial(pathway_hsa, name_gene):
    gl.logger.debug('Get information about the input gene "%s" in the pathway "%s"' % (name_gene, pathway_hsa))

    filename = os.path.join(os.getcwd(), 'database', 'pathways', 'xml', pathway_hsa + '.xml.gz')

    with gzip.open(filename, "rb") as f:
        content = f.read().decode('utf-8')
        mydoc = minidom.parseString(content)

        # GENES
        entry = mydoc.getElementsByTagName('entry')

        for elem in entry:
            string_check = name_gene + ','
            if elem.attributes['type'].value == 'gene' and string_check in \
                    elem.getElementsByTagName('graphics')[0].attributes['name'].value:

                gl.logger.debug('Find gene information "%s"! hsa: "%s" and url: "%s"' % (
                    name_gene, elem.attributes['name'].value, elem.attributes['link'].value))

                return elem.attributes['name'].value, elem.attributes['link'].value

        print(gl.COLORS['red'] + gl.COLORS['bold'] + "The gene entered not exit into pathway selected! " +
              gl.COLORS['end_line'])
        gl.logger.error('The gene "%s" entered not exit into pathway selected "%s"' % (name_gene, pathway_hsa))

        exit(1)


def download_file(url, pathfile, filename):
    if not os.path.exists(os.path.join(pathfile, filename)):
        #gl.logger.warning('The file "%s" does not exist in the selected path ("%s")' % (filename, pathfile))
        try:
            req = requests.get(url)

            with gzip.open(os.path.join(pathfile, filename), "wb") as f:
                f.write(req.content)
            #gl.logger.warning('The file "%s" was created in the selected path ("%s")' % (filename, pathfile))
        except requests.exceptions.ConnectionError:
            print("Connection refused from KEGG")
            #gl.logger.error('Connection refused from KEGG')
            exit(1)
    #else:
    #    gl.logger.warning('The file "%s" exists in the selected path ("%s")' % (filename, pathfile))


def download_read_html(url):
    filename = url.split('?')[1].replace('+', '_').replace(':', '')
    filename = hashlib.md5(filename.encode("utf-8")).hexdigest() + '.html.gz'

    # scarico file html con tutti i pathway collegati a geni passati nell'url
    download_file(url, os.path.join(os.getcwd(), 'database', 'pathways', 'html'), filename)

    # open pathway of the genes in html
    with gzip.open(os.path.join(os.getcwd(), 'database', 'pathways', 'html', filename), "rb") as f:
        gl.logger.debug('Reading of the html file obtained, downloaded from the url "%s"' % url)

        content = f.read().decode('utf-8')

        soup = bs.BeautifulSoup(content, 'html.parser')

        list_pathway = list()
        for link in soup.findAll('a'):
            a_tag = str(link.get('href'))
            if "/kegg-bin/show_pathway" in a_tag:
                a_tag = a_tag.split("+")
                a_tag = a_tag[0].split("?")
                list_pathway.append(a_tag[1])

        list_pathway = list(set(list_pathway))  # generate a unique elements

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


def read_kgml(hop, pathway_hsa, name_gene_start, hsa_gene_start, path, occu):
    # print('[%s][hop: %s][pathway: %s][gene: %s]\n' % (i, hop, pathway_hsa, name_gene_start))
    #gl.logger.debug('[hop: %s] Reading the gene "%s" contained in the pathway "%s"' % (
     #   hop, name_gene_start, pathway_hsa))

    filename = os.path.join(os.getcwd(), 'database', 'pathways', 'xml', pathway_hsa + '.xml.gz')

    with gzip.open(filename, "rb") as f:

        content = f.read().decode('utf-8')
        mydoc = minidom.parseString(content)

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

        # print('[%s][hop: %s][pathway: %s][gene: %s]list_genes_this_pathway: %s\n' % (i, hop, pathway_hsa,
        # name_gene_start, list_genes_this_pathway))

        # cerco gli id all'interno della mappa per quello specifico hsa(gene)
        # perchè in ogni pathway è diverso, anche se è lo stesso gene
        # se li trova restituisce tuple di diversi id che fanno riferimento al gene
        list_ids_gene_input = search_id_to_hsa(list_genes_this_pathway, hsa_gene_start)

        # print('[%s][hop: %s][pathway: %s][gene: %s]list_ids_gene_input: %s\n' % (i, hop, pathway_hsa,
        # name_gene_start, list_ids_gene_input))

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

                        # print('[%s][hop: %s][pathway: %s][gene: %s]list_gene_relation: %s\n' % (
                        # i, hop, pathway_hsa, name_gene_start, list_gene_relation))

                        # ONLY DEBUG
                        # if len(list_gene_relation) == 0:
                        #    print('mi blocco')

                        # print('FIND! @ ID: ' + elem.attributes['entry1'].value + ' - name: ' + list_gene_relation[
                        # 0][2])

                        # verifico se esiste una relazione effettivamente
                        if len(list_gene_relation) > 0:
                            #gl.logger.debug('[hop: %s] Found a direct relationship between "%s" and "%s" in the '
                            #                'pathway "%s" with a relationship "%s"' % (
                            #                    hop, name_gene_start, list_gene_relation[0][2], pathway_hsa,
                            #                    elem.attributes['type'].value))

                            # in concat_multiple_subtype, concateno tutti i subtype della relazione analizzata
                            row = {
                                'hop': hop,
                                'name_start': name_gene_start,
                                'hsa_start': hsa_gene_start,
                                'name_end': list_gene_relation[0][2],
                                'hsa_end': list_gene_relation[0][1],
                                'url_gene_end': list_gene_relation[0][3],
                                'relation': elem.attributes['type'].value,
                                'type_rel': concat_multiple_subtype(elem.getElementsByTagName('subtype')),
                                'pathway_origin': pathway_hsa,
                                'path': path + '/' + list_gene_relation[0][2],
                                'occurrences': occu
                            }
                            list_rows.append(row)
        return list_rows


def analysis_hop_n(hop, gene, gene_hsa, pathway_this_gene, path, occu):
    # print('[i: %s][hop: %s][gene: %s][hsa: %s][pathway: %s]\n' % (iteration, hop, gene, gene_hsa, pathway_this_gene))

    # download_xml(pathway_this_gene)
    download_file('http://rest.kegg.jp/get/' + pathway_this_gene + '/kgml',
                  os.path.join(os.getcwd(), 'database', 'pathways', 'xml'),
                  pathway_this_gene + '.xml.gz')

    list_rows = read_kgml(hop, pathway_this_gene, gene, gene_hsa, path, occu)

    return list_rows


def unified(list_rows_returned):
    gl.logger.info('Add the results obtained by the threads in parallel in the main data frame')
    for row in list_rows_returned:
        if row is not None:
            for cell in row:
                gl.DF_TREE = gl.DF_TREE.append(cell, ignore_index=True)


def get_info_row_duplicated(i, hop, df_filtered, gene):
    grouped_df = (df_filtered[df_filtered['name_end'] == gene]).groupby("name_start")

    list_to_do_df = list()
    iter = 0
    for key, group in grouped_df:
        # key = gene
        # group = group of this gene
        # print('[i_par: %s][key: %s][iter_for: %s] group2: %s\n' % (i, key, iter, group.iloc[0]['path']))
        # print('[i_par: %s][key: %s][iter_for: %s] shape_row: %s\n' % (i, key, iter, group.shape[0]))

        hsa_end_refactor = ' '.join(group["hsa_end"].tolist())
        hsa_end_refactor = sorted(set(hsa_end_refactor.split(' ')))
        hsa_end_refactor = ' '.join(hsa_end_refactor)
        # print('[i_par: %s][key: %s][iter_for: %s] hsa_end_refactor2: %s\n' % (i, key, iter, hsa_end_refactor))

        url_gene_end_refactor = 'https://www.kegg.jp/dbget-bin/www_bget?' + hsa_end_refactor.replace(' ', '+')
        occurences_calculated = group.iloc[0]['occurrences'] * group.shape[0]

        # print('[i_par: %s][key: %s][iter_for: %s] url_gene_end_refactor: %s\n' % (i, key, iter, url_gene_end_refactor))

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
            '§§'.join(group['pathway_origin'].tolist()),
            occurences_calculated
        )
        )
        iter = iter + 1

    return list_to_do_df


def clean_update_row_duplicates(list_to_do_df):
    for row in list_to_do_df:
        for cell in row:
            # print(row, '\n')

            # aggiorno la riga selezionata, con le nuove stringhe nelle 4 colonne
            cols = ['hsa_end', 'url_gene_end', 'relation', 'type_rel', 'pathway_origin', 'occurrences']

            gl.DF_TREE.loc[cell[0], cols] = [cell[2], cell[3], cell[4], cell[5], cell[6], cell[7]]

            # rimuovo le righe selezionate
            if len(cell[1]) > 0:
                gl.DF_TREE = gl.DF_TREE.drop(cell[1])
