import globals as gl
import time
from utility import read_config, clear_previous_results, load_last_csv, create_zip, check_database
from draw import draw_json_run
from analysis import run_analysis
import os

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

print("----- CHECK DATABASE -----")
check_database()

print("----- START ANALYSIS -----")
run_analysis(starting_depth)
print("----- END ANALYSIS -----")

print("----- START GENERATE OUTPUT -----")
draw_json_run(os.path.join(os.getcwd(), 'export_data', 'df_resulted.csv'))
print("----- END GENERATE OUTPUT -----")

print("----- START GENERATE ZIPFILE -----")
create_zip()
print("----- END GENERATE ZIPFILE -----")

m, s = divmod(time.time() - start_time, 60)
print(f"----- DONE EXECUTION ({round(m)} mins, {round(s)} secs) -----")
