import pandas as pd
import requests
import io
import os
import gzip
import multiprocessing as mlp
from joblib import Parallel, delayed


def download_gene_from_KEGG(i, hsa):
    url = "https://www.kegg.jp/dbget-bin/www_bget?%s" % hsa
    filename = "%s.html" % hsa
    if not os.path.exists(os.path.join(os.getcwd(), 'database', 'pathways', 'html', filename)):
        try:
            req = requests.get(url)

            with gzip.open(os.path.join(os.getcwd(), 'database', 'pathways', 'html', filename), "wb") as f:
                f.write(req.content)

        except requests.exceptions.ConnectionError:
            print("ERROR: Connection refused from KEGG")
            exit(1)
    print("Downloaded the file \"%s\" with index: %s" % (hsa, i))


def get_list_genes_from_KEGG():
    url = "http://rest.kegg.jp/list/hsa"
    try:
        req = requests.get(url).content.decode('utf-8')
        df_temp = pd.read_csv(io.StringIO(req), sep="	", names=['hsa', 'gene'])

        return df_temp['hsa']
    except requests.exceptions.ConnectionError:
        print("ERROR: Connection refused from KEGG for get list of the hsa")
        exit(1)


def download_pathway_from_KEGG(i, hsa):
    url = "http://rest.kegg.jp/get/%s/kgml" % hsa
    filename = "%s.xml.gz" % hsa
    if not os.path.exists(os.path.join(os.getcwd(), 'database', 'pathways', 'xml', filename)):
        try:
            req = requests.get(url)

            with gzip.open(os.path.join(os.getcwd(), 'database', 'pathways', 'xml', filename), "wb") as f:
                f.write(req.content)

        except requests.exceptions.ConnectionError:
            print("ERROR: Connection refused from KEGG")
            exit(1)

        print("Downloaded the file \"pathway %s\" with index: %s" % (hsa, i))


def get_list_pathways_from_KEGG():
    url = "http://rest.kegg.jp/list/pathway"
    try:
        req = requests.get(url).content.decode('utf-8')
        df_temp = pd.read_csv(io.StringIO(req), sep="	", names=['map', 'name'])

        df_temp['map'] = df_temp.map.str.replace('path:map', 'hsa')

        return df_temp['map']
    except requests.exceptions.ConnectionError:
        print("ERROR: Connection refused from KEGG for get list of the pathways")
        exit(1)


#df_genes_hsa = get_list_genes_from_KEGG()

#Parallel(n_jobs=mlp.cpu_count())(delayed(download_gene_from_KEGG)(i, hsa) for i, hsa in df_genes_hsa.iteritems())


df_pathways_hsa = get_list_pathways_from_KEGG()
Parallel(n_jobs=mlp.cpu_count())(delayed(download_pathway_from_KEGG)(i, hsa) for i, hsa in df_pathways_hsa.iteritems())
