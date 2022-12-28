import argparse
import time
from multiprocessing import cpu_count
from os.path import join as pathjoin, exists as pathexists
from os import getcwd, makedirs
from shutil import rmtree


def header():
    t = '========================================================\n' \
        '=          PETAL – ParallEl paThways AnaLyzer          =\n' \
        '=                        v2.0                          =\n' \
        '=          Last update:     2021/06/29                 =\n' \
        '=          database update: 2021/06/20                 =\n' \
        '========================================================\n' \
        '=          E-mail: giuseppe.sgroi@unict.it             =\n' \
        '========================================================\n' \
        '=          PETAL is licensed under CC BY-NC-SA 4.0     =\n' \
        '========================================================'
    print(t)


def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {'false', 'f', '0', 'no', 'n'}:
        return False
    elif value.lower() in {'true', 't', '1', 'yes', 'y'}:
        return True
    raise ValueError(f'{value} is not a valid boolean value')


def read_args():
    parser = argparse.ArgumentParser(description='PETAL requires the setting of five mandatory '
                                                 'and optional input parameters')
    subprasers = parser.add_subparsers(dest='command')

    menu1 = subprasers.add_parser('analysis', help='Analysis: the in-depth analysis will be launched '
                                                   'with the generation of the textual and graphic tree')

    menu1.add_argument('--load', type=str_to_bool, nargs='?', const=True, default=False)

    menu1.add_argument('-p', '--pathway_id', type=str, help='Biological pathway id (in hsa format) '
                                                            '– Default value = hsa04010', default='hsa04010')
    menu1.add_argument('-g', '--gene_id', type=str, help='Starting gene (id) present in the selected pathway '
                                                         '– Default value = hsa:5594 (MAPK1)', default='hsa:5594')
    menu1.add_argument('-d', '--depth', type=int, help='Maximum search depth of the analysis – Default value = 2',
                       default=2)
    menu1.add_argument('-c', '--cpu', type=int, help='(optional) Maximum number of CPUs used '
                                                     'during the analysis – Default value = 0',
                       choices=range(0, cpu_count()), default=2)

    menu2 = subprasers.add_parser('filter', help='Filter: it allows to filter through a target gene the '
                                                 'results previously obtained with the generation of a textual tree')
    menu2.add_argument('-t', '--targets', type=str, help='Target gene used to filter the data', required=True)
    menu2.add_argument('-c', '--cpu', type=int, help='(optional) Maximum number of CPUs used '
                                                     'during the analysis – Default value = 1',
                       choices=range(0, cpu_count()), default=1)

    args = parser.parse_args()

    if hasattr(args, 'cpu') and args.cpu == 0:
        args.cpu = cpu_count()

    print(f'\t\u279C {args}')

    return args


def clear_results():
    print('\t\u26D4\u0020\u0020Removing files in progress')
    path = pathjoin(getcwd(), 'export_data')

    if pathexists(path):
        rmtree(path)

    makedirs(path)

    time.sleep(2)

    print('\t\u2705\u0020\u0020DONE')
