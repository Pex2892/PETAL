import globals as gl
from utility import download_file, download_read_html, set_progress_bar, export_data_for_level
import os
import gzip
from joblib import Parallel, delayed
from xml.dom.minidom import parseString


def run_analysis():
    for deep in range(1, gl.deep_input + 1):
        if deep == 1:
            # download initial pathway
            download_file('http://rest.kegg.jp/get/' + gl.pathway_input + '/kgml',
                              os.path.join(os.getcwd(), 'database', 'pathways', 'xml'), gl.pathway_input + '.xml.gz')

            # get info first gene from gene name
            hsa_finded, url_finded = get_info_gene_initial(gl.pathway_input, gl.gene_input)

            # set globals variables
            gl.gene_input_hsa = hsa_finded
            gl.gene_input_url = url_finded

            # read initial pathway, create and add genes to csv
            list_rows_df_returned = read_kgml(deep, gl.pathway_input, gl.gene_input, hsa_finded, gl.gene_input, 1)

            # add n genes found to the dataframe
            unified([list_rows_df_returned])

            # retrive other list pathways in reference to initial pathway
            list_pathways_this_gene = download_read_html(url_finded)

            # The pathway set as input from the config file is removed
            if gl.pathway_input in list_pathways_this_gene:
                list_pathways_this_gene.remove(gl.pathway_input)

            # process single gene on each CPUs available
            list_rows_df_returned = Parallel(n_jobs=gl.num_cores_input)(delayed(analysis_deep_n)(
                deep, gl.gene_input, hsa_finded, pathway_this_gene, gl.gene_input, 1)
                                                                        for pathway_this_gene in set_progress_bar(
                '[Deep: %d]' % deep, str(len(list_pathways_this_gene)))(list_pathways_this_gene))

            unified(list_rows_df_returned)

        else:
            # Retrieve the genes found at depth-1, avoiding the input gene
            df_genes_resulted = (
                gl.DF_TREE[(gl.DF_TREE['deep'] == deep - 1) & (gl.DF_TREE['name_son'] != gl.gene_input)])

            for index, row in set_progress_bar(
                    '[Deep: %d]' % deep, str(df_genes_resulted.shape[0]))(df_genes_resulted.iterrows()):
                # Return a list of pathways about the gene passed in input
                list_pathways_this_gene = download_read_html(row['url_kegg_son'])

                # The pathway set as input from the config file is removed, so as to avoid an endless loop
                if gl.pathway_input in list_pathways_this_gene:
                    list_pathways_this_gene.remove(gl.pathway_input)

                # process single gene on each CPUs available
                list_rows_df_returned = Parallel(n_jobs=gl.num_cores_input)(
                    delayed(analysis_deep_n)(
                        deep, row['name_son'], row['hsa_son'], pathway_this_gene,
                        row['fullpath'], row['occurrences']) for pathway_this_gene in list_pathways_this_gene
                )

                unified(list_rows_df_returned)

        # ----- DROP DUPLICATES -----

        # Duplicates of the same level are extracted and sorted in alphabetical order
        df_genes_this_level = (gl.DF_TREE[gl.DF_TREE['deep'] == deep])
        df_duplicated_filtered = df_genes_this_level[df_genes_this_level.duplicated(
            subset=['name_son'], keep=False)].sort_values('name_son')

        # The names of the genes that are duplicated are recovered
        list_name_genes_duplicated = df_duplicated_filtered.name_son.unique()

        # process single gene on each CPUs available
        list_rows_to_do_df_returned = Parallel(n_jobs=gl.num_cores_input)(
            delayed(get_info_row_duplicated)(df_duplicated_filtered, gene_duplicate)
            for gene_duplicate in list_name_genes_duplicated
        )

        # The number of occurrences of the found links is updated and the duplicates will be deleted
        clean_update_row_duplicates(list_rows_to_do_df_returned)

        # Row indexes are reset, because they are no longer sequential due to the elimination of duplicates
        gl.DF_TREE = gl.DF_TREE.reset_index(drop=True)

        export_data_for_level(deep)

        # ----- DROP DUPLICATES -----


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


def search_id_to_hsa(list_genes_this_pathway, hsa_gene):
    return [item for item in list_genes_this_pathway if hsa_gene == item[1]]


def search_gene_to_id(list_genes_this_pathway, id_gene):
    return [item for item in list_genes_this_pathway if id_gene in item]


def concat_multiple_subtype(list_subtype):
    # All the "subtypes" of the selected relationship are concatenated
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
                # Genes are saved are of the "gene" type
                list_genes_this_pathway.append((elem.attributes['id'].value, elem.attributes['name'].value,
                                                elem2.attributes['name'].value.split(',')[0],
                                                elem.attributes['link'].value))

        # All the ids contained in the pathway are recovered, because a gene can have
        # different ids in different pathways or in the same pathway
        list_ids_gene_input = search_id_to_hsa(list_genes_this_pathway, hsa_gene_start)

        list_rows = list()
        if len(list_ids_gene_input) > 0:

            # You scroll through the relationships found within the pathway
            for elem in relation:
                # Scroll the list of ids found
                for id_gene in list_ids_gene_input:
                    # Check that "entry2" has at least one match in the id list
                    if elem.attributes['entry2'].value == id_gene[0]:

                        # We try to find a correspondence between gene and your id. If it returns zero,
                        # it indicates that that gene does not exist or is group or compund type and not gene type
                        list_gene_relation = search_gene_to_id(list_genes_this_pathway, elem.attributes['entry1'].value)

                        # Checks if at least one connection has been found
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

        # row to update (take the first one by updating all fields)
        # list of rows to be removed, keeping only one
        # the hsa_end of the gene is created
        # url_end of the gene is created
        # concatenation of all relations based on the source pathways
        # concatenation of all type_rel based on the source pathway
        # concatenation of all path_origin
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
            # The selected row is updated, in the 4 specified columns, with the new information
            cols = ['hsa_son', 'url_kegg_son', 'relation', 'type_rel', 'pathway_of_origin', 'occurrences']

            gl.DF_TREE.loc[cell[0], cols] = [cell[2], cell[3], cell[4], cell[5], cell[6], cell[7]]

            # The selected rows are removed, because they are duplicated
            if len(cell[1]) > 0:
                gl.DF_TREE = gl.DF_TREE.drop(cell[1])
