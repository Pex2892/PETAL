import time
from utility import header, read_args, clear_results, create_zip
from kegg import KEGG
from analysis import Analysis
from draw import Draw
from filter import Filter


start_time = time.time()

header()

print(">>> READ PARAMETERS <<<")
args = read_args()

print(">>> CHECK DATABASE <<<")
obj = KEGG('https://github.com/Pex2892/PETAL/releases/download/v1.2/only_database.zip', 172800)
print(">>> LOADED LIST OF HUMAN GENES <<<")

if args.command == 'analysis':
    obj = Analysis(args, obj)

    if args.load:
        print(">>> LOADING PREVIOUS RESULTS <<<")
        obj.load_results()
    else:
        print(">>> DELETING PREVIOUS RESULTS <<<")
        clear_results()

    print(">>> STARTING – ANALYSIS <<<")
    obj.run()
    r = obj.__repr__()
    print(">>> ENDING – ANALYSIS <<<")

    print(">>> STARTING – GENERATION OF THE TREE <<<")
    obj = Draw(args)
    obj.from_analysis(r)
    print(">>> ENDING – GENERATION OF THE TREE <<<")

    print(">>> STARTING – GENERATION OF THE ZIP FILE <<<")
    create_zip(f'analysis_{args.pathway}_{args.gene}_{args.depth}')
    print(">>> ENDING – GENERATION OF THE ZIP FILE <<<")

elif args.command == 'filter':
    print(">>> STARTING – FILTER <<<")
    obj = Filter(args, obj)
    r = obj.run()

    print(">>> STARTING – GENERATION OF THE TEXT TREE <<<")
    obj = Draw(args)
    obj.from_filter(r)
    print(">>> ENDING – GENERATION OF THE TEXT TREE <<<")

    print(">>> STARTING – GENERATION OF THE ZIP FILE <<<")
    create_zip('analysis_with_filters')
    print(">>> ENDING – GENERATION OF THE ZIP FILE <<<")

    print(">>> ENDING – FILTER <<<")


m, s = divmod(time.time() - start_time, 60)
print(f"----- DONE EXECUTION ({round(m)} mins, {round(s)} secs) -----")
