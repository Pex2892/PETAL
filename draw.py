import os
import shutil
import pandas as pd
from joblib import Parallel, delayed
from anytree import Node, RenderTree, ContRoundStyle
from anytree.exporter import JsonExporter
from shutil import copyfile, copytree
from utility import pb

COLS = ['depth', 's_gene', 's_gene_hsa', 'e_gene', 'e_gene_hsa',
        'e_gene_url', 'isoforms', 'relation', 'type_rel',
        'pathway_origin', 'fullpath', 'occurrences']


class Draw:
    def __init__(self, args):
        self.input_cpu = args.cpu

    def from_analysis(self, info_gene: list):
        df = pd.read_csv(os.path.join(os.getcwd(), 'export_data', 'df_resulted.csv'), sep=";", names=COLS)

        path_lists = df['fullpath'].str.split('/').values.tolist()

        info_lists = df[['depth', 'e_gene_hsa', 'e_gene_url', 'isoforms',
                         'relation', 'type_rel', 'pathway_origin', 'occurrences']].values.tolist()

        tree = self.list_to_anytree(path_lists, info_lists, info_gene)

        path_tree = os.path.join(os.getcwd(), 'export_data', 'demo_radialtree')
        if os.path.exists(path_tree):
            shutil.rmtree(path_tree)
        os.makedirs(path_tree)

        self.export_tree_in_json(tree, path_tree)

        self.export_tree_in_txt(tree, path_tree)

        copyfile(os.path.join(os.getcwd(), 'static', 'html', 'index.html'), os.path.join(path_tree, 'index.html'))
        copyfile(os.path.join(os.getcwd(), 'static', 'html', 'help.html'), os.path.join(path_tree, 'help.html'))
        copytree(os.path.join(os.getcwd(), 'static', 'assets'), os.path.join(path_tree, 'assets'))

    def from_filter(self, p_list: list):
        Parallel(n_jobs=self.input_cpu, backend='threading') \
            (delayed(self.from_filter_par)(p) for p in pb('[Draw filters]', len(p_list))(p_list))

    def from_filter_par(self, p: str):
        df = pd.read_csv(os.path.join(p, 'df_filtered.csv'), sep=";", names=COLS)
        path_lists = df['fullpath'].str.split('/').values.tolist()
        info_lists = df[['depth', 'e_gene_hsa', 'e_gene_url', 'isoforms',
                         'relation', 'type_rel', 'pathway_origin', 'occurrences']].values.tolist()

        tree = self.list_to_anytree(path_lists, info_lists)

        self.export_tree_in_txt(tree, p)

    def list_to_anytree(self, lst, lst2, info_root=None):
        root_name = lst[0][0]

        if info_root is None:
            root_node = Node(name=root_name)
        else:
            root_node = Node(name=root_name, hsa=info_root[0], url=info_root[2], info='NaN',
                             occurrences=0, deep=0, isoforms='NaN')

        for branch, i in zip(lst, lst2):
            parent_node = root_node
            assert branch[0] == parent_node.name

            for cur_node_name in branch[1:]:
                cur_node = next(
                    (node for node in parent_node.children if node.name == cur_node_name),
                    None,
                )
                if cur_node is None:
                    cur_node = Node(name=cur_node_name, hsa=i[1], url=i[2], info=self.concat_info(i[4], i[5], i[6]),
                                    occurrences=i[7], deep=i[0], isoforms=str(i[3]), parent=parent_node)
                parent_node = cur_node
        return root_node

    def concat_info(self, rel, type_rel, patwhay):
        # AGGIUNGERE ISOFORME e USARE i nel for e non lo zip
        rel_arr = rel.split('§§')
        type_rel_arr = type_rel.split('§§')
        patwhay_arr = patwhay.split('§§')

        str_info = ' - '.join([f'{c} | {a} | {b}' for a, b, c in zip(rel_arr, type_rel_arr, patwhay_arr)])
        return str_info

    def export_tree_in_json(self, tree, path):
        f = open(os.path.join(path, 'data-flare.json'), 'w')
        exporter = JsonExporter(indent=4)
        exporter.write(tree, f)

    def export_tree_in_txt(self, tree, path):
        tree_txt = ''
        for pre, _, node in RenderTree(tree, style=ContRoundStyle()):
            tree_txt += f'{pre}{node.name}\n'

        with open(os.path.join(path, 'tree.txt'), 'w') as f:
            f.write(tree_txt)
