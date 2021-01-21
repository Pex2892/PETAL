import argparse
import datetime
import os
import subprocess as sb
import time
from draw import draw_from_filter
from utility import read_list_homo_sapiens_genes, check_gene_and_alias, get_gene_info_from_name, create_zip


print("----- START FILTER -----")
start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--gtarget', help='Target gene used to filter the data', required=True)
args = parser.parse_args()

target_lists = args.gtarget.split(',')
# print(target_lists)

filepath = os.path.join(os.getcwd(), 'export_data', 'df_resulted.csv')

for target in target_lists:
    list_all_genes = read_list_homo_sapiens_genes()

    gene_info_target = get_gene_info_from_name(target, list_all_genes)
    # print(gene_info_target)

    if gene_info_target is not None:
        gene_target = check_gene_and_alias(target, gene_info_target[1])

        filterdir_path = os.path.join(os.getcwd(), 'export_data',
                                      f'filter_{target}_{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}')
        os.makedirs(filterdir_path)
        filename = os.path.join(filterdir_path, 'df_filtered.csv')

        try:
            sb.check_output(f'grep -wE \'{gene_target}|{gene_info_target[0]}\' {filepath} > {filename}', shell=True)
        except sb.CalledProcessError as e:
            print('The filter did not find any results!')
        else:
            draw_from_filter(filterdir_path)
    else:
        print('The inserted target gene does not exist! Try to verify manually.')

print("----- END FILTER -----")

print("----- START GENERATE ZIPFILE -----")
create_zip(f'analysis_with_filters')
print("----- END GENERATE ZIPFILE -----")

m, s = divmod(time.time() - start_time, 60)
print(f"----- DONE EXECUTION ({round(m)} mins, {round(s)} secs) -----")