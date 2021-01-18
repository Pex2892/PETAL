import os
from utility import read_list_homo_sapiens_genes, check_gene_and_alias, get_gene_info_from_name
import subprocess as sb

gene_input = 'MVCD1'

filepath = os.path.join(os.getcwd(), 'export_data', 'df_resulted.csv')
filename = os.path.join(os.getcwd(), 'export_data', f'df_filtered.csv')

list_all_genes = read_list_homo_sapiens_genes()

gene_info = get_gene_info_from_name(gene_input, list_all_genes)
print(gene_info)

gene_input = check_gene_and_alias(gene_input, gene_info[1])

try:
    sb.check_output(f'grep -wE \'{gene_input}|{gene_info[0]}\' {filepath} > {filename}', shell=True)
except sb.CalledProcessError as e:
    print('Nessun valore è stato trovato!')
    exit()


