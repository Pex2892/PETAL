from pandas import read_csv
from os.path import join as pathjoin, exists as pathexists
from os import getcwd, makedirs
from anytree import Node, RenderTree, ContRoundStyle
from shutil import copyfile, copytree, rmtree
from anytree.exporter import JsonExporter


class Tree:
    def __init__(self, args, kegg):
        self.args = args
        self.kegg = kegg

    def run(self, flag):
        if flag:
            self._from_analysis()

    def _from_analysis(self):
        df = read_csv(pathjoin(getcwd(), 'export_data', 'analysis_all_depths.csv'), sep='\t')

        root = [self.args.gene_id, self.kegg.df_genes.at[self.args.gene_id, 'name'].split(',', 1)[0],
                self.kegg.df_genes.at[self.args.gene_id, 'url']]
        path_lists = df['fullpath'].str.split('/').values.tolist()

        info_lists = df[['depth', 'ending_gene_id', 'ending_gene_name', 'ending_isoforms_id',
                         'ending_isoforms_name', 'relation', 'subtype', 'reference_pathway', 'occurrences']].values.tolist()

        res_tree = self._list_to_anytree(path_lists, info_lists, root)

        path_tree = pathjoin(getcwd(), 'export_data', 'demo_radialtree')
        if pathexists(path_tree):
            rmtree(path_tree)
        makedirs(path_tree)

        self._export_tree_in_json(res_tree, path_tree)

        self._export_tree_in_txt(res_tree, path_tree)

        copyfile(pathjoin(getcwd(), 'static', 'html', 'index.html'), pathjoin(path_tree, 'index.html'))
        copyfile(pathjoin(getcwd(), 'static', 'html', 'help.html'), pathjoin(path_tree, 'help.html'))
        copytree(pathjoin(getcwd(), 'static', 'assets'), pathjoin(path_tree, 'assets'))

        # print()

    def _list_to_anytree(self, lst, lst2, info_root=None):
        root_name = lst[0][0]

        if info_root is None:
            root_node = Node(name=root_name)
        else:
            root_node = Node(name=root_name, hsa=info_root[0], url=info_root[2], info='NaN',
                             occurrences=0, deep=0, isoforms='NaN')

        # print(lst)
        for branch, i in zip(lst, lst2):
            # print(i)
            parent_node = root_node
            assert branch[0] == parent_node.name

            for cur_node_name in branch[1:]:
                cur_node = next((node for node in parent_node.children if node.name == cur_node_name), None,)
                if cur_node is None:
                    cur_node = Node(name=cur_node_name, hsa=i[1],
                                    url=f'https://www.genome.jp/dbget-bin/www_bget?{i[1]}',
                                    info=f'{i[6]} ({i[5]})',
                                    occurrences=i[8], deep=i[0], isoforms=f'{i[3]}', parent=parent_node)
                parent_node = cur_node
        return root_node

    def _export_tree_in_json(self, tree, path):
        f = open(pathjoin(path, 'data-flare.json'), 'w')
        exporter = JsonExporter(indent=4)
        exporter.write(tree, f)

    def _export_tree_in_txt(self, tree, path):
        tree_txt = ''
        for pre, _, node in RenderTree(tree, style=ContRoundStyle()):
            # print(node)
            tree_txt += f'{pre}{node.name}\n'

        with open(pathjoin(path, 'tree.txt'), 'w') as f:
            f.write(tree_txt)
