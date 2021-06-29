import time
from os.path import join as pathjoin, exists as pathexists
from os import getcwd

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from itertools import chain

class Analysis:
    def __init__(self, args, kegg):
        self.args = args
        self.kegg = kegg
        self.df_results_prev = None
        self.starting_depth = 1

        if self.args.load is False:
            self._checks()

    def run(self):
        self._breadth_analysis()

    def _checks(self):
        time.sleep(0)
        self._check_exist_pathway()

        self._check_exist_gene()

        self._check_exist_gene_into_pathway()

    def _check_exist_pathway(self):
        print('\t\u26D4\u0020\u0020Checking that the pathway id exists in the database')
        if self.args.pathway_id not in self.kegg.df_pathways.index:
            print(f'\t\u279C The pathway-id "{self.args.pathway_id}" was not found. Please try to verify the parameter.')
            exit()
        time.sleep(0)
        print('\t\u2705\u0020\u0020DONE')

    def _check_exist_gene(self):
        print('\t\u26D4\u0020\u0020Checking that the gene id exists in the database')
        if self.args.gene_id not in self.kegg.df_genes.index:
            print(f'\t\u279C The gene-id "{self.args.gene_id}" was not found. Please try to verify the parameter.')
            exit()
        time.sleep(0)
        print('\t\u2705\u0020\u0020DONE')

    def _check_exist_gene_into_pathway(self):
        print('\t\u26D4\u0020\u0020Checks that the gene is present within the chosen pathway')
        if self.args.gene_id not in self.kegg.df_pathways.at[self.args.pathway_id, 'kgml']:
            print(f'\t\u279C The gene-id "{self.args.gene_id}" was not found within the "{self.args.pathway_id}" pathway.')
            exit()
        time.sleep(0)
        print('\t\u2705\u0020\u0020DONE')

    def _breadth_analysis(self):

        if self.starting_depth == 1:

            fullpath = self.kegg.df_genes.at[self.args.gene_id, 'name'].split(',', 1)[0]

            results = self.kegg.read_kgml(self.args.pathway_id, self.args.gene_id, fullpath)
            df_results = pd.DataFrame(results)
            df_results['depth'] = self.starting_depth
        else:
            subdf = self.df_results_prev[self.df_results_prev['ending_gene_id'] != self.args.gene_id]

            results = []
            # devo linearizzare questo for ed eliminare il parallelo dal for esterno
            # troppo lento così come è adesso
            '''for i in subdf.index:
                gene_id = self.df_results_prev.at[i, 'ending_gene_id']
                list_pathways = self.kegg.df_genes.at[gene_id, 'pathways'].split(',')

                if self.args.pathway_id in list_pathways:
                    list_pathways.remove(self.args.pathway_id)

                more_info = [{'pathway_id': list_pathways[j],
                              'fullpath': self.df_results_prev.at[i, 'fullpath']} for j in range(0, len(list_pathways))]

                r = Parallel(n_jobs=self.args.cpu, verbose=5)(delayed(self.kegg.read_kgml)(
                    list_pathways[j], gene_id, more_info[j]) for j in range(0, len(list_pathways)))
                flat_list = [item for sublist in r for item in sublist]
                results.append(flat_list)'''

            to_elab = []
            for i in subdf.index:
                gene_id = self.df_results_prev.at[i, 'ending_gene_id']
                list_pathways = self.kegg.df_genes.at[gene_id, 'pathways'].split(',')

                if self.args.pathway_id in list_pathways:
                    list_pathways.remove(self.args.pathway_id)

                for j in range(0, len(list_pathways)):
                    to_elab.append({
                        'gene_id': gene_id,
                        'pathway_id': list_pathways[j],
                        'fullpath': self.df_results_prev.at[i, 'fullpath'],
                    })

            r = Parallel(n_jobs=self.args.cpu, verbose=5, backend='threading')(delayed(self.kegg.read_kgml)(
                to_elab[i]['pathway_id'], to_elab[i]['gene_id'], to_elab[i]['fullpath'])
                                                          for i in range(0, len(to_elab)))
            flat_list = list(chain(*r))

            df_results = pd.DataFrame(flat_list)
            df_results['depth'] = self.starting_depth

        self._duplicate_results(df_results)

        self._export()

        print(f'\t\u2705\u0020\u0020Depth {self.starting_depth}: completed')

        self.starting_depth += 1
        if self.starting_depth <= self.args.depth:
            self.run()

    def _duplicate_results(self, df_results):
        duplicate_rows = df_results[df_results.duplicated(df_results.columns.tolist()[0:-1], keep=False)]

        if duplicate_rows.shape[0] > 0:
            temp_df = duplicate_rows.groupby(df_results.columns.tolist()[0:-1])

            for name, group in temp_df:
                df_results.at[group.index[0], 'occurrences'] *= group.shape[0]
                df_results.loc[group.index[1:], :] = np.nan

        df_results.dropna(inplace=True)

        self.df_results_prev = df_results

    def _export(self):
        self.df_results_prev.to_csv(pathjoin(getcwd(), 'export_data', f'analysis_depth_{self.starting_depth}.csv'),
                                    sep='\t', header=True, index=False)

        self.df_results_prev.to_csv(pathjoin(getcwd(), 'export_data', f'analysis_all_depths.csv'), sep='\t', mode='a',
                                    header=not pathexists(pathjoin(getcwd(), 'export_data', f'analysis_all_depths.csv')),
                                    index=False)
