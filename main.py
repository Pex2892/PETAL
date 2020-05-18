import time
from utility import read_config, clear_previous_results, check_pathway_update_history, export_data
from draw import draw_json_run
from analysis import run_analysis

# --------------- INITIAL START TIME --------------
start_time = time.time()

# -------------- INITIAL MAIN --------------
print("----- INITIAL SHELL PARAMETERS -----")
read_config()

print("----- CLEAN PREVIOUS RESULTS -----")
clear_previous_results()

print("----- CHECK UPDATED PATHWAYS -----")
check_pathway_update_history('https://www.genome.jp/kegg/docs/upd_map.html')

print("----- START ANALYSIS -----")
run_analysis()

print("----- EXPORT DATA -----")
export_data()

print("----- START GENERATE OUTPUT -----")

draw_json_run()

m, s = divmod(time.time() - start_time, 60)
print("----- DONE EXECUTION (%d mins, %d secs) -----" % (round(m), round(s)))
