import argparse
import os
import shutil
import requests
import subprocess as sb
import multiprocessing as mlp
import gzip
from progressbar import Bar, Counter, ETA, Percentage, ProgressBar


def header():
    t = '========================================================\n' \
         '=          PETAL – ParallEl paThways AnaLyzer          =\n' \
         '=                        v1.3                          =\n' \
         '=          Last update:     2021/02/10                 =\n' \
         '=          database update: 2020/12/24                 =\n' \
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

    analysis = subprasers.add_parser('analysis', help='Analysis: the in-depth analysis will be launched '
                                                      'with the generation of the textual and graphic tree')

    analysis.add_argument('--load', type=str_to_bool, nargs='?', const=True, default=False)

    analysis.add_argument('-p', '--pathway', type=str, help='Biological pathway (in hsa format) '
                                                            '– Default value = hsa04010', default='hsa04010')
    analysis.add_argument('-g', '--gene', type=str, help='Starting gene present in the selected pathway '
                                                         '– Default value = MAPK1', default='MAPK1')
    analysis.add_argument('-d', '--depth', type=int, help='Maximum search depth of the analysis '
                                                          '– Default value = 2', default=2)
    analysis.add_argument('-c', '--cpu', type=int, help='(optional) Maximum number of CPUs used '
                                                        'during the analysis – Default value = 0',
                          choices=range(0, mlp.cpu_count()), default=0)

    filter = subprasers.add_parser('filter', help='Filter: it allows to filter through a target gene the '
                                                  'results previously obtained with the generation of a textual tree')
    filter.add_argument('-t', '--targets', type=str, help='Target gene used to filter the data', required=True)
    filter.add_argument('-c', '--cpu', type=int, help='(optional) Maximum number of CPUs used '
                                                      'during the analysis – Default value = 0',
                        choices=range(0, mlp.cpu_count()), default=0)

    args = parser.parse_args()

    if args.cpu == 0:
        args.cpu = mlp.cpu_count()

    print(args)

    return args


def clear_results():
    path = os.path.join(os.getcwd(), 'export_data')

    if os.path.exists(path):
        shutil.rmtree(path)

    os.makedirs(path)


def download_file(url, filename, destination, message='ERROR: Connection refused'):
    try:
        req = requests.get(url)

        if req.text:
            with gzip.open(os.path.join(destination, filename), "wb") as f:
                f.write(req.content)
            return True
        else:
            print(f"The url \"{url}\" returned empty content!")
            return False
    except requests.exceptions.ConnectionError:
        print(message)
        exit(1)


def pb(action: str, max_elem: int):
    pb = ProgressBar(widgets=[action, ' ', Percentage(), ' (', Counter(), ' of ',
                              str(max_elem), ') ', Bar('#'), ' ', ETA()], maxval=max_elem)
    return pb


def create_zip(f_name: str):
    root_path = os.getcwd()
    path = os.path.join(os.getcwd(), 'export_data')

    try:
        sb.check_output(f'cd {path} && '
                        f'zip -r {f_name}.zip . -x \'*.DS_Store\' -x \'*.zip\' > /dev/null && '
                        f'cd {root_path}', shell=True)
    except sb.CalledProcessError:
        print('An error occurred while compressing the results.')
