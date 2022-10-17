import globals as gl
import time
import os
from analysis import run_analysis
from utility import read_params, clear_previous_results, load_last_csv, create_zip, check_database, header
from draw import draw_from_analysis


# --------------- INITIAL START TIME --------------
start_time = time.time()

header()

# -------------- INITIAL MAIN --------------
print("----- INITIAL SHELL PARAMETERS -----")
read_params()

if gl.mode_input == 0:
    print("----- CLEAN PREVIOUS RESULTS -----")
    clear_previous_results()
    starting_depth = 1
else:
    print("----- LOAD LAST RESULTS (CSV) SAVED -----")
    starting_depth = load_last_csv()

print("----- CHECK DATABASE -----")
check_database()

print("----- START ANALYSIS -----")
run_analysis(starting_depth)
print("----- END ANALYSIS -----")

print("----- START GENERATE OUTPUT -----")
draw_from_analysis(
    [gl.gene_input_hsa, gl.gene_input, gl.gene_input_url],
    os.path.join(os.getcwd(), 'export_data')
)
print("----- END GENERATE OUTPUT -----")

print("----- START GENERATE ZIPFILE -----")
create_zip(f'analysis_{gl.pathway_input}_{gl.gene_input}_{gl.deep_input}')
print("----- END GENERATE ZIPFILE -----")

m, s = divmod(time.time() - start_time, 60)
print(f"---- DONE EXECUTION ({round(m)} mins, {round(s)} secs) -----")
