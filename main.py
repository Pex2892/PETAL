import globals as gl
import time
import subprocess
from joblib import Parallel, delayed
import utility as utl
from generator import *

# --------------- INITIAL START TIME --------------
start_time = time.time()

# -------------- INITIAL MAIN --------------
# read initial parameters
utl.read_config()

# delete all results of previous execution
utl.clear_previous_results()

# download and return list of the update pathway
utl.check_pathway_update_history('https://www.genome.jp/kegg/docs/upd_map.html')

# hop+1 = because it has to analyze the last level
for level_actual in range(1, gl.hop_input + 1):
    print(gl.COLORS['pink'] + gl.COLORS['bold'] + "--- START LEVEL %s ---" % level_actual + gl.COLORS['end_line'])

    if level_actual == 1:
        # download initial pathway
        utl.download_file('http://rest.kegg.jp/get/' + gl.pathway_input + '/kgml',
                          os.path.join(os.getcwd(), 'database', 'pathways', 'xml'),
                          gl.pathway_input + '.xml.gz')

        # get info first gene from hsa name of pathway
        hsa_gene_input_finded, url_pathway_gene_input_finded = utl.get_info_gene_initial(gl.pathway_input,
                                                                                         gl.gene_input)

        # read initial pathway, create and add genes to csv
        list_rows_df_returned = utl.read_kgml(level_actual, gl.pathway_input, gl.gene_input, hsa_gene_input_finded)

        # unifico i primi n geni che sono direttamente connessi
        utl.unified([list_rows_df_returned])

        # retrive other list pathways in reference to initial pathway
        list_pathways_this_gene = utl.download_read_html(url_pathway_gene_input_finded)

        # rimuovo il pathway di origine, cioè passato in input dal file di config
        if gl.pathway_input in list_pathways_this_gene:
            list_pathways_this_gene.remove(gl.pathway_input)

        # process single gene on each CPUs available
        list_rows_df_returned = Parallel(n_jobs=gl.num_cores_input, backend='threading')(delayed(utl.analysis_hop_n)(
            level_actual, gl.gene_input, hsa_gene_input_finded,
            pathway_this_gene) for pathway_this_gene in list_pathways_this_gene)

        utl.unified(list_rows_df_returned)

    else:
        df_genes_resulted = (gl.DF_TREE[gl.DF_TREE['hop'] == level_actual - 1])

        for index, row in df_genes_resulted.iterrows():
            # print('[%s] SONO ARRIVATOOOOOO: %s ' % (level_actual, index))

            # ottengo la lista di pathway in riferimento al gene che sto passando
            list_pathways_this_gene = utl.download_read_html(row['url_gene_end'])

            # rimuovo il pathway di origine, cioè passato in input dal file di config così evito un loop continuo
            if gl.pathway_input in list_pathways_this_gene:
                list_pathways_this_gene.remove(gl.pathway_input)

            # process single gene on each CPUs available
            list_rows_df_returned = Parallel(n_jobs=gl.num_cores_input, backend='threading')(
                delayed(utl.analysis_hop_n)(
                    level_actual, row['name_end'], row['hsa_end'],
                    pathway_this_gene) for pathway_this_gene in list_pathways_this_gene)

            utl.unified(list_rows_df_returned)

    # ----- DROP DUPLICATES -----

    # estraggo i duplicati dello stesso livello e ordinati in ordine alfabetico
    df_genes_this_level = (gl.DF_TREE[gl.DF_TREE['hop'] == level_actual])
    df_duplicated_filtered = df_genes_this_level[df_genes_this_level.duplicated(
        subset=['name_end'], keep=False)].sort_values('name_end')

    # lista con i nomi dei geni duplicati
    list_name_genes_duplicated = df_duplicated_filtered.name_end.unique()

    # process single gene on each CPUs available
    list_rows_to_do_df_returned = Parallel(n_jobs=gl.num_cores_input, backend='threading')(
        delayed(utl.get_info_row_duplicated)(
            df_duplicated_filtered, gene_duplicate) for gene_duplicate in list_name_genes_duplicated)

    # aggiorno e elimino le righe del dataframe
    utl.clean_update_row_duplicates(list_rows_to_do_df_returned)

    # resetto l'indice di riga, perchè non più sequenziali dovuto alle eliminazioni delle righe
    gl.DF_TREE = gl.DF_TREE.reset_index(drop=True)

    # ----- DROP DUPLICATES -----

    print(gl.COLORS['pink'] + gl.COLORS['bold'] + "--- END LEVEL %s ---" % level_actual + gl.COLORS['end_line'])

print(gl.DF_TREE)

# subprocess.call('sh concatenate.sh', shell=True)
gl.DF_TREE.to_csv('results/results_level.csv', sep=';', header=False, index=False)


# genero il file di output file output_text.txt
# init_generate_output(gl.hop)

# genero il file di output in json dal file output_text.txt
# output_json()

m, s = divmod(time.time() - start_time, 60)

print(gl.COLORS['red'] + gl.COLORS['bold'] + "--- %s minutes and %s seconds ---" % (round(m), round(s)) + gl.COLORS[
    'end_line'])
