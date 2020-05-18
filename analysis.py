import globals as gl
import utility as utl
import os
from joblib import Parallel, delayed


def run_analysis():
    for deep in range(1, gl.deep_input + 1):
        if deep == 1:
            # download initial pathway
            utl.download_file('http://rest.kegg.jp/get/' + gl.pathway_input + '/kgml',
                              os.path.join(os.getcwd(), 'database', 'pathways', 'xml'), gl.pathway_input + '.xml.gz')

            # get info first gene from gene name
            hsa_finded, url_finded = utl.get_info_gene_initial(gl.pathway_input, gl.gene_input)

            # set globals variables
            gl.gene_input_hsa = hsa_finded
            gl.gene_input_url = url_finded

            # read initial pathway, create and add genes to csv
            list_rows_df_returned = utl.read_kgml(deep, gl.pathway_input, gl.gene_input, hsa_finded, gl.gene_input, 1)

            # add n genes found to the dataframe
            utl.unified([list_rows_df_returned])

            # retrive other list pathways in reference to initial pathway
            list_pathways_this_gene = utl.download_read_html(url_finded)

            # rimuovo il pathway di origine, cioè passato in input dal file di config
            if gl.pathway_input in list_pathways_this_gene:
                list_pathways_this_gene.remove(gl.pathway_input)

            # process single gene on each CPUs available
            list_rows_df_returned = Parallel(n_jobs=gl.num_cores_input)(delayed(utl.analysis_deep_n)(
                deep, gl.gene_input, hsa_finded, pathway_this_gene, gl.gene_input, 1)
                                                                        for pathway_this_gene in utl.set_progress_bar(
                '[Deep: %d]' % deep, str(len(list_pathways_this_gene)))(list_pathways_this_gene))

            utl.unified(list_rows_df_returned)

        else:
            # Retrieve the genes found at depth-1, avoiding the input gene
            df_genes_resulted = (
            gl.DF_TREE[(gl.DF_TREE['deep'] == deep - 1) & (gl.DF_TREE['name_son'] != gl.gene_input)])

            for index, row in utl.set_progress_bar('[Deep: %d]' % deep, str(df_genes_resulted.shape[0])) \
                        (df_genes_resulted.iterrows()):
                # ottengo la lista di pathway in riferimento al gene che sto passando
                list_pathways_this_gene = utl.download_read_html(row['url_kegg_son'])

                # rimuovo il pathway di origine, cioè passato in input dal file di config così evito un loop continuo
                if gl.pathway_input in list_pathways_this_gene:
                    list_pathways_this_gene.remove(gl.pathway_input)

                # process single gene on each CPUs available
                list_rows_df_returned = Parallel(n_jobs=gl.num_cores_input)(
                    delayed(utl.analysis_deep_n)(deep, row['name_son'], row['hsa_son'], pathway_this_gene,
                                                 row['fullpath'], row['occurrences']) for pathway_this_gene in
                    list_pathways_this_gene
                )

                utl.unified(list_rows_df_returned)

        # ----- DROP DUPLICATES -----

        # estraggo i duplicati dello stesso livello e ordinati in ordine alfabetico
        df_genes_this_level = (gl.DF_TREE[gl.DF_TREE['deep'] == deep])
        df_duplicated_filtered = df_genes_this_level[df_genes_this_level.duplicated(
            subset=['name_son'], keep=False)].sort_values('name_son')

        # lista con i nomi dei geni duplicati
        list_name_genes_duplicated = df_duplicated_filtered.name_son.unique()

        # process single gene on each CPUs available
        list_rows_to_do_df_returned = Parallel(n_jobs=gl.num_cores_input)(
            delayed(utl.get_info_row_duplicated)(df_duplicated_filtered, gene_duplicate)
            for gene_duplicate in list_name_genes_duplicated
        )

        # aggiorno e elimino le righe del dataframe
        utl.clean_update_row_duplicates(list_rows_to_do_df_returned)

        # resetto l'indice di riga, perchè non più sequenziali dovuto alle eliminazioni delle righe
        gl.DF_TREE = gl.DF_TREE.reset_index(drop=True)

        # ----- DROP DUPLICATES -----