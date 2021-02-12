import os
import re
import glob
import pandas as pd
from joblib import Parallel, delayed
from utility import pb


class Analysis:
    def __init__(self, args, obj_kegg):
        self.input_pathway = args.pathway
        self.input_gene = args.gene
        self.input_depth = args.depth
        self.input_cpu = args.cpu
        self.KEGG = obj_kegg
        self.input_gene_hsa = None
        self.input_gene_url = None
        self.starting_depth = None
        self.results_per_depth = None

        self.preanalysis()

    def __repr__(self):
        return repr([self.input_gene_hsa, self.input_gene, self.input_gene_url])

    def preanalysis(self):
        t = self.KEGG.get_gene_from_name(self.input_gene, True)

        self.input_gene = self.KEGG.check_alias(self.input_gene, t)
        self.input_gene_hsa = t[0]
        self.input_gene_url = t[2]

        self.KEGG.check_exist_gene_in_pathway(self.input_pathway, self.input_gene)

        self.starting_depth = 1
        cols = ['depth', 's_gene', 's_gene_hsa', 'e_gene', 'e_gene_hsa',
                'e_gene_url', 'isoforms', 'relation', 'type_rel',
                'pathway_origin', 'fullpath', 'occurrences']
        self.results_per_depth = pd.DataFrame(columns=cols)

    def run(self):
        if self.starting_depth == 1:
            # retrive other list pathways in reference to initial pathway
            pathways_list = self.KEGG.read_gene_txt(self.input_gene_hsa)

            # print(pathways_list, len(pathways_list))
            r = Parallel(n_jobs=self.input_cpu, backend='threading') \
                (delayed(self.KEGG.read_kgml)(
                    [self.starting_depth, pathway, self.input_gene, self.input_gene_hsa, self.input_gene, 1])
                 for pathway in pb(f'[Depth: {self.starting_depth}]', len(pathways_list))(pathways_list))

            Parallel(n_jobs=self.input_cpu, backend='sequential')(delayed(self.append_results)(row) for row in r)
            del r
        else:
            # Retrieve the genes found at depth-1, avoiding the input gene
            df_filter = self.results_per_depth[(self.results_per_depth['depth'] == self.starting_depth - 1)
                                               & self.results_per_depth['e_gene'] != self.input_gene]

            for index, row in pb(f'[Depth: {self.starting_depth}]', df_filter.shape[0])(df_filter.iterrows()):
                # Return a list of pathways about the gene passed in input
                pathways_list = self.KEGG.read_gene_txt(row['e_gene_hsa'])

                # process single gene on each CPUs available
                r = Parallel(n_jobs=self.input_cpu, backend='threading')(delayed(self.KEGG.read_kgml)(
                    [self.starting_depth, p, row['e_gene'], row['e_gene_hsa'], row['fullpath'], row['occurrences']])
                    for p in pathways_list)

                Parallel(n_jobs=self.input_cpu, backend='sequential')(delayed(self.append_results)(row) for row in r)
                del r

        self.preprocessing_results()

        print(f'Depth {self.starting_depth}: Completed.')

        self.results_per_depth = self.results_per_depth[(self.results_per_depth['depth'] == self.starting_depth)]

        self.export_results_per_depth(self.starting_depth)
        print(f'Depth {self.starting_depth}: CSV saved.')

        print('---------------')

        self.starting_depth += 1
        if self.starting_depth <= self.input_depth:
            self.run()

    def append_results(self, row):
        # Add the results obtained by the threads in parallel in the main data frame
        # for row in list_rows_returned:
        if row is not None:
            self.results_per_depth = self.results_per_depth.append(row, ignore_index=True)

    def preprocessing_results(self):
        df = self.results_per_depth[self.results_per_depth.duplicated(
            subset=['depth', 's_gene', 's_gene_hsa', 'e_gene',
                    'e_gene_hsa', 'e_gene_url'], keep=False)].sort_values('e_gene')
        # print(self.results_per_depth_filtered)

        # The names of the genes that are duplicated are recovered
        name_genes_list = df.e_gene.unique()
        # print(name_genes_list)

        # process single gene on each CPUs available
        r = Parallel(n_jobs=self.input_cpu)(delayed(self.get_info_row_dupl)(df, g) for g in name_genes_list)

        # The number of occurrences of the found links is updated and the duplicates will be deleted
        self.clean_update_row_duplicates(r)

        # print(self.results_per_depth)

        # Row indexes are reset, because they are no longer sequential due to the elimination of duplicates
        self.results_per_depth = self.results_per_depth.reset_index(drop=True)

        # print(self.results_per_depth)

    def get_info_row_dupl(self, df, gene):
        df_g = (df[df['e_gene'] == gene]).groupby('s_gene')

        list_to_do_df = list()
        for key, group in df_g:
            occurences_calculated = group.iloc[0]['occurrences'] * group.shape[0]

            # row to update (take the first one by updating all fields)
            # list of rows to be removed, keeping only one
            # concatenation of all isoforms based on the source pathways
            # concatenation of all relations based on the source pathways
            # concatenation of all type_rel based on the source pathway
            # concatenation of all path_origin
            list_to_do_df.append(
                (
                    group.index[0],
                    list(filter(group.index[0].__ne__, group.index.values.tolist())),
                    '§§'.join(group['isoforms'].tolist()),
                    '§§'.join(group['relation'].tolist()),
                    '§§'.join(group['type_rel'].tolist()),
                    '§§'.join(group['pathway_origin'].tolist()),
                    occurences_calculated
                )
            )
        return list_to_do_df

    def clean_update_row_duplicates(self, list_to_do_df):
        for row in list_to_do_df:
            for cell in row:
                # The selected row is updated, in the 4 specified columns, with the new information
                cols = ['isoforms', 'relation', 'type_rel', 'pathway_origin', 'occurrences']

                self.results_per_depth.loc[cell[0], cols] = [cell[2], cell[3], cell[4], cell[5], cell[6]]

                # The selected rows are removed, because they are duplicated
                if len(cell[1]) > 0:
                    self.results_per_depth = self.results_per_depth.drop(cell[1])

    def export_results_per_depth(self, depth):
        self.results_per_depth.to_csv(os.path.join(os.getcwd(), 'export_data', f'df_resulted_depth_{depth}.csv'),
                                      sep=';', header=False, index=False)

        self.results_per_depth.to_csv(os.path.join(os.getcwd(), 'export_data', f'df_resulted.csv'), sep=';', mode='a',
                                      header=False, index=False)

    def load_results(self):
        list_files = glob.glob(os.path.join(os.getcwd(), 'export_data', 'df_resulted_depth_*.csv'))
        list_files.sort(key=lambda f: int(re.sub(r'\D', '', f)))  # path sorted

        # path, f_name, depth
        info_last_csv = (list_files[-1], os.path.basename(list_files[-1]),
                         int(re.compile(r'\d+').findall(os.path.basename(list_files[-1]))[0]))

        # print(info_last_csv, self.input_depth)

        if info_last_csv[2] >= self.input_depth:
            print('>>>>>> The selected maximum depth has already been analyzed.'
                  'You will find the results in the "exporta_data" folder.')
            exit(1)
        elif info_last_csv[2] < self.input_depth:
            self.results_per_depth = pd.concat(
                [self.results_per_depth, pd.read_csv(info_last_csv[0], sep=";", names=self.results_per_depth.columns)])
            self.starting_depth = info_last_csv[2] + 1

            print(f'>>>>>> The analysis will re-start from a depth of {self.starting_depth}, '
                  f'because previous analyses have already been performed (i.e. depths up to {self.starting_depth-1}).')
