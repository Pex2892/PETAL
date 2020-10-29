import globals as gl
import time
from utility import read_config, clear_previous_results, check_pathway_update_history, load_last_csv, create_zip, API_KEGG_get_list_human_genes
from draw import draw_json_run
from analysis import run_analysis

# --------------- INITIAL START TIME --------------
start_time = time.time()

# -------------- INITIAL MAIN --------------
print("----- INITIAL SHELL PARAMETERS -----")
read_config()

if gl.mode_input == 0:
    print("----- CLEAN PREVIOUS RESULTS -----")
    clear_previous_results()
    starting_depth = 1
else:
    print("----- LOAD LAST RESULTS (CSV) SAVED -----")
    starting_depth = load_last_csv()

print("----- CHECK UPDATED PATHWAYS -----")
check_pathway_update_history('https://www.genome.jp/kegg/docs/upd_map.html')

print("-----  LOAD LIST OF HUMAN GENES -----")
API_KEGG_get_list_human_genes()

print("----- START ANALYSIS -----")
run_analysis(starting_depth)
print("----- END ANALYSIS -----")

print("----- START GENERATE OUTPUT -----")
draw_json_run()
print("----- END GENERATE OUTPUT -----")

print("----- START GENERATE ZIPFILE -----")
create_zip()
print("----- END GENERATE ZIPFILE -----")

m, s = divmod(time.time() - start_time, 60)
print("----- DONE EXECUTION (%d mins, %d secs) -----" % (round(m), round(s)))
