#!/usr/bin/python3
# -*- coding: utf-8 -*-

if __name__ == '__main__':
    import os
    import argparse
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    from cpartition import x2wp, CProfiles

    rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': 'Arial',
                     'font.size': 13,
                     'mathtext.fontset': 'stix'})

    y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
             Si=5.02504411E-2, Fe=9.4414085022e-1)

    parser = argparse.ArgumentParser()
    parser.add_argument('basenames', nargs='+')
    parser.add_argument('-S', '--save', action='store_true',
                        help='Save plot')
    args = parser.parse_args()

    for basename in args.basenames:
        cprofiles = CProfiles(basename)

        print(cprofiles.basename)

        ax = cprofiles.plot_colormap(mirror=True,
                                     func=lambda x: x2wp(x, y=y),
                                     vmin=0, vmax=1.8)
        ax.set_xlabel(u'Position (μm)')
        ax.set_ylabel('Time (s)')
        ax.set_title(cprofiles.basename)

        if args.save:
            plt.savefig(os.path.join('img', cprofiles.basename + '.png'), dpi=150)
            plt.close()

    plt.show()
