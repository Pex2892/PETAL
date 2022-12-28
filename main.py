import time
import os
import http.server
from utility import header, read_args, clear_results
from database import Kegg
from analysis import Analysis
from tree import Tree

start_time = time.time()

header()

print(f'\u276F\u276F\u276F READ ARGUMENTS \u276E\u276E\u276E')
args = read_args()

print('\u276F\u276F\u276F CHECKING THE KEGG DATABASE \u276E\u276E\u276E')
kegg = Kegg()
kegg.run()

if args.command == 'analysis' or args.load is False:
    print('\u276F\u276F\u276F DELETING PREVIOUS RESULTS \u276E\u276E\u276E')
    clear_results()
    print('\u276F\u276F\u276F CHECKING INPUT PARAMETERS \u276E\u276E\u276E')
    analysis = Analysis(args, kegg)
    print('\u276F\u276F\u276F ANALYSIS \u276E\u276E\u276E')
    analysis.run()

    tree = Tree(args, kegg)
    tree.run(True)


m, s = divmod(time.time() - start_time, 60)
print(f"----- DONE EXECUTION ({round(m)} mins, {round(s)} secs) -----")

os.chdir('export_data/demo_radialtree')
http.server.test(
    HandlerClass=http.server.SimpleHTTPRequestHandler,
    bind='127.0.0.1',
    port=8080)

