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

                            print('cancello il pathway vecchio (locale), ma aggiornato su kegg: %s' %
                                  ('hsa' + merged_pathway))
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


def download_read_html(url):
    # download and read of the html page containing all the genes, passed as parameters in the url

    filename = url.split('?')[1].replace('+', '_').replace(':', '')
    filename = hashlib.md5(filename.encode("utf-8")).hexdigest() + '.html.gz'

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


def set_progress_bar(action, max_elem):
    pb = ProgressBar(widgets=[action, ' ', Percentage(), ' (', Counter(), ' of ',
                              max_elem, ') ', Bar('#'), ' ', ETA()], maxval=int(max_elem))
    return pb


def export_data():
    gl.DF_TREE.to_csv(os.path.join(os.getcwd(), 'export_data', 'df_resulted.csv'), sep=';', header=False, index=False)
