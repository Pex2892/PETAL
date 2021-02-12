import os
import datetime
import subprocess as sb


class Filter:
    def __init__(self, args, obj_kegg):
        self.targets = args.targets.split(',')
        self.input_cpu = args.cpu
        self.KEGG = obj_kegg

    def run(self):
        f_path = os.path.join(os.getcwd(), 'export_data', 'df_resulted.csv')

        r_draw = []
        for t in self.targets:
            info = self.KEGG.get_gene_from_name(t)

            # check exist gene
            if info is not None:
                t = self.KEGG.check_alias(t, info)

                filterdir = os.path.join(os.getcwd(), 'export_data',
                                         f'filter_{t}_{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}')
                os.makedirs(filterdir)
                f_name = os.path.join(filterdir, 'df_filtered.csv')

                # print(t, info[0], f_path, f_name)
                try:
                    sb.check_output(f'grep -wE \'{t}|{info[0]}\' {f_path} > {f_name}',
                                    shell=True)
                except sb.CalledProcessError:
                    os.remove(f_name)
                    message = f'The filter "{t}" did not find any results.'
                    with open(os.path.join(filterdir, f'{t}_README.txt'), 'w') as f:
                        f.write(message)
                    print('>>>>>>', message)
                else:
                    print(f'>>>>>> The filter "{t}" was completed successfully.')
                    r_draw.append(filterdir)
        return r_draw
