import time
import multiprocessing
import subprocess
from joblib import Parallel, delayed
from utility import *
from generator import *
from Bcolors import Bcolors as Bc


# --------------- INITIAL START TIME --------------
start_time = time.time()


# -------------- INITIAL CONFIG --------------
pd.set_option('display.max_colwidth', -1) # mi serve a non troncare le info salvate all'interno del DF
col = ['hop', 'name_start', 'hsa_start', 'name_end', 'hsa_end', 'url_hsa_end', 'relation', 'type_rel', 'pathway_origin']
df_tree = pd.DataFrame(columns=col)


# -------------- INITIAL SHELL PARAMETERS --------------
pathway_hsa, name_gene, hop = command_line()

print(Bc.WARNING + Bc.BOLD + "--- INITIAL SHELL PARAMETERS ---" + Bc.ENDC)
print(Bc.OKBLUE + Bc.BOLD + "Pathway: %s" % pathway_hsa + Bc.ENDC)
print(Bc.OKBLUE + Bc.BOLD + "Gene: %s" % name_gene + Bc.ENDC)
print(Bc.OKBLUE + Bc.BOLD + "Hop: %s" % hop + Bc.ENDC)


hsa_gene = 'hsa:5594'  # non deve essere statico, ma prelevarlo da questo url:
                       # https://www.genome.jp/dbget-bin/www_bget?pathway+hsa04010
# https://www.genome.jp/dbget-bin/www_bget?hsa:5594+hsa:5595



# -------------- INITIAL MAIN --------------

# delete all results of previous execution
clear_results()

# set the number of the CPUs
num_cores = multiprocessing.cpu_count()

# hop+1 = because it has to analyze the last level
for level_actual in range(1, hop+1):
    print(Bc.HEADER + Bc.BOLD + "--- START LEVEL %s ---" % level_actual + Bc.ENDC)

    if level_actual == 1:
        download_xml(pathway_hsa)
        read_kgml(level_actual, pathway_hsa, name_gene, hsa_gene)

        list_pathways_gene_actual_first_level = [x for x in download_read_html('https://www.genome.jp/dbget-bin/www_bget?hsa:5594+hsa:5595') if x != pathway_hsa]

        inputs = range(0, len(list_pathways_gene_actual_first_level))

        Parallel(n_jobs=num_cores)(delayed(execute_i)('https://www.genome.jp/dbget-bin/www_bget?hsa:5594+hsa:5595', pathway_hsa, level_actual, name_gene, hsa_gene, list_pathways_gene_actual_first_level[i]) for i in inputs)
    else:
        df_genes_level_prev = (df_tree[df_tree['hop'] == level_actual - 1])

        inputs = range(0, df_genes_level_prev.shape[0])

        Parallel(n_jobs=num_cores)(delayed(execute_i)(df_genes_level_prev.iloc[i,7], pathway_hsa, level_actual, df_genes_level_prev.iloc[i,3], df_genes_level_prev.iloc[i,4], None) for i in inputs)

    # carico il csv del livello attuale con i duplicati
    df_tree = pd.read_csv('results/results_level'+str(level_actual)+'.csv', sep=';', header=None, names=col)

    # elimino i duplicati dal livello attuale e ri-salvo il file
    df_tree = df_tree.drop_duplicates(subset='name_end')
    df_tree.to_csv('results/results_level'+str(level_actual)+'.csv', sep=';', header=False, index=False)

    print(Bc.HEADER + Bc.BOLD + "--- END LEVEL %s ---" % level_actual + Bc.ENDC)


subprocess.call('sh concatenate.sh', shell=True)

# genero il file di output file output_text.txt
init_generate_output(hop)

# genero il file di output in json dal file output_text.txt
output_json()

m, s = divmod(time.time() - start_time, 60)

print(Bc.FAIL + Bc.BOLD + "--- %s minutes and %s seconds ---" % (round(m), round(s)) + Bc.ENDC)
