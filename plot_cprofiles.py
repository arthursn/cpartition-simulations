#!/usr/bin/python3
# -*- coding: utf-8 -*-

if __name__ == '__main__':
    import sys
    import os
    import numpy as np
    from fnmatch import fnmatch
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    from cpartition import x2wp, CProfiles
    import argparse

    rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': 'Arial',
                     'font.size': 13})

    y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
             Si=5.02504411E-2, Fe=9.4414085022e-1)

    dictloc = dict(ur='upper right', ul='upper left',
                   ll='lower left', lr='lower right',
                   best='best')

    parser = argparse.ArgumentParser()
    parser.add_argument('basenames', nargs='*')
    parser.add_argument('-s', '--show', action='store_true')
    parser.add_argument('-S', '--save', action='store_true')

    parser.add_argument('-d', '--dir',
                        default='/home/arthur/tese/img/cpartition/cprofiles')
    parser.add_argument('-e', '--ext', default='.svg')
    parser.add_argument('-a', '--append', default='')
    parser.add_argument('-f', '--figsize', type=float, nargs=2, default=[6, 4])

    parser.add_argument('-m', '--mirror', action='store_true')
    parser.add_argument('-x', '--xlim', type=float, nargs=2, default=[None, None])
    parser.add_argument('-y', '--ylim', type=float, nargs=2, default=[None, None])

    parser.add_argument('-t', '--time', type=float, nargs='*', required=True)
    parser.add_argument('-l', '--label', action='store_true')
    parser.add_argument('-T', '--tracking', action='store_true')
    parser.add_argument('-L', '--loc', default='best')

    args = parser.parse_args()

    for basename in args.basenames:
        if 'mart' in basename:
            labels = [('aus1', r'$\gamma_1$', 1),
                      ('aus2', r'$\gamma_2$', 1),
                      ('aust', r'$\gamma$', 1),
                      ('mart', r"$\alpha'$", 1)]
            tracked_interfaces = []
        else:
            labels = [('aus1', r'$\gamma_1$', -1),
                      ('aus2', r'$\gamma_2$', -1),
                      ('aust', r'$\gamma$', -1),
                      ('mart', r"$\alpha'$", -1),
                      ('fer1', r'$\alpha_{b1}$', 1),
                      ('fer2', r'$\alpha_{b2}$', 1)]
            tracked_interfaces = [('aus1.sn', 'aus1.cin'),
                                  ('aus2.s0', 'aus2.ci0'),
                                  ('aus2.sn', 'aus2.cin')]

        try:
            fig, ax = plt.subplots(figsize=args.figsize)
            cprofiles = CProfiles(basename)
            cprofiles.plot_cprofiles(ax=ax, tlist=args.time,
                                     mirror=args.mirror,
                                     func=lambda x: x2wp(x, y=y))

        except:
            print('Failed to plot "{}"'.format(basename))
            plt.close()
        else:
            ax.set_xlim(*args.xlim)
            ax.set_ylim(*args.ylim)
            ax.set_xlabel(u'Posição (μm)')
            ax.set_ylabel('Teor de carbono (%)')
            ax.legend(loc=dictloc[args.loc], fancybox=False)

            if args.tracking:
                cprofiles.plot_locus_interface(tracked_interfaces, ax=ax,
                                               mirror=args.mirror,
                                               color='k', ls='--', lw=.8,
                                               func=lambda x: x2wp(x, y=y), label='')

            if args.label:
                cprofiles.label_phases(ax=ax, t=args.time[0],
                                       labels=labels,
                                       mirror=args.mirror, size=12)

            if args.save:
                fname = os.path.join(args.dir,
                                     cprofiles.basename + args.append + args.ext)
                plt.savefig(fname, bbox_inches='tight')
                os.system('svg2pdf ' + fname)

    if args.show:
        plt.show()
    else:
        plt.close('all')
