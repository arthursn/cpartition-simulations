#!/usr/bin/python3

# -*- coding: utf-8 -*-

if __name__ == '__main__':
    import sys
    import os
    import numpy as np
    import pandas as pd
    import matplotlib
    import matplotlib.pyplot as plt
    from cpartition import FCC, BCC, Interface, WBs, CProfiles, x2wp
    from itertools import cycle
    import argparse

    colorcycle = cycle(
        matplotlib.rcParams['axes.prop_cycle'].by_key()['color'])

    y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
             Si=5.02504411E-2, Fe=9.4414085022e-1)

    # For calculating WBs
    T_C = 375.
    aust = FCC(T_C=T_C, tdata='thermo/FoFo/TCFE8/375-fcc.txt')
    ferr = BCC(T_C=T_C, tdata='thermo/FoFo/TCFE8/375-bcc.txt', E=WBs(T_C))
    intf = Interface(domain1=aust, domain2=ferr, type_int='mobile.eq')
    # WBs composition
    _, cwbs = intf.comp()
    cwbs = x2wp(cwbs, y=y)

    parser = argparse.ArgumentParser()
    parser.add_argument('basenames', nargs='+')
    parser.add_argument('-s', '--silent', action='store_false')
    parser.add_argument('-S', '--save', action='store_true')

    parser.add_argument('-d', '--dir',
                        default='/home/arthur/tese/img/cpartition/cprofiles')
    parser.add_argument('-e', '--ext', default='.svg')
    parser.add_argument('-a', '--append', default='')
    parser.add_argument('-f', '--figsize', type=float, nargs=2, default=[6, 4])

    parser.add_argument('-m', '--mirror', action='store_true')
    parser.add_argument('-x', '--xlim', type=float,
                        nargs=2, default=[None, None])
    parser.add_argument('-y', '--ylim', type=float,
                        nargs=2, default=[None, None])

    parser.add_argument('-t', '--time', type=float, nargs='+', required=True)
    parser.add_argument('-l', '--label', action='store_true')
    parser.add_argument('-T', '--tracking', action='store_true')
    parser.add_argument('-L', '--loc', default='best')

    parser.add_argument('--special', action='store_true')
    parser.add_argument('--all', action='store_true')
    parser.add_argument('--wbs', action='store_true')

    args = parser.parse_args()

    if args.all:
        args.time.append(args.time.copy())

    ncol = 2
    nrow = int(np.ceil(len(args.time)/ncol))

    for basename in args.basenames:
        if 'mart' in basename:
            labels = [('aus1', r'$\gamma_1$', 1),
                      ('aus2', r'$\gamma_2$', 1),
                      ('aust', r'$\gamma$', 1)]

            if 'CCEpara' in basename or 'CCEortho' in basename or 'mu' in basename:
                labels += [('mart', r"$\alpha' + \theta$", 1)]
            else:
                labels += [('mart', r"$\alpha'$", 1)]
        else:
            labels = [('aus1', r'$\gamma_1$', -1),
                      ('aus2', r'$\gamma_2$', -1),
                      ('aust', r'$\gamma$', -1),
                      ('fer1', r'$\alpha_{b1}$', 1),
                      ('fer2', r'$\alpha_{b2}$', 1)]

            if 'CCEpara' in basename or 'CCEortho' in basename or 'mu' in basename:
                labels += [('mart', r"$\alpha' + \theta$", -1)]
            else:
                labels += [('mart', r"$\alpha'$", -1)]

        cprofiles = CProfiles(basename, 'C_profiles')

        if args.tracking:
            pos = pd.read_table(
                'pos_extremities/{}.txt'.format(cprofiles.basename), sep=' ')
            ci = pd.read_table(
                'C_extremities/{}.txt'.format(cprofiles.basename), sep=' ')

        fig, axes = plt.subplots(nrow, ncol, figsize=(4*ncol, 3*nrow))
        # plt.subplots_adjust(wspace=.5, hspace=.5)
        axes = axes.ravel()

        for i, (t, ax) in enumerate(zip(args.time, axes)):
            kw = dict(lw=1)
            if not isinstance(t, list):
                t = [t]
                kw['color'] = next(colorcycle)

            cprofiles.plot_cprofiles(ax=ax, mirror=True,
                                     func=lambda x: x2wp(x, y=y),
                                     tlist=t, **kw)
            if i == 0:
                j = cprofiles.where_tlist(t, [])
                if len(j) > 0:
                    j = j[0]
                    idx, = np.where(cprofiles.ss[j] == 'aus1')
                    idx = idx[0]
                    zmax = 2*cprofiles.zz[j][-1] - cprofiles.zz[j][idx]
                    cmax = cprofiles.cc[j][idx]

            if args.tracking:
                ax.plot(pos['aus1.sn'], x2wp(ci['aus1.cin'], y=y), 'k:')

            ax.set_ylim(*args.ylim)
            ax.set_xlim(*args.xlim)

            if len(t) == 1:
                ax.text(.01, .85, s='{:g} s'.format(t[0]),
                        weight='bold', transform=ax.transAxes)

                if args.wbs:
                    ax.axhline(cwbs, color='k', ls='--', lw=.8)
                    ax.text(1.1, cwbs + .1, s='WBs', ha='right')

                # special axis dimensions
                if args.special:
                    if i > 1:
                        ax.set_ylim(-.1, 2.5)
                    else:
                        ax.set_ylim(-.1)

                cprofiles.label_phases(ax, t, labels=labels, mirror=True)
            else:
                ax.text(.01, .85, s='Todos',
                        weight='bold', transform=ax.transAxes)

            if i == len(axes) - ncol:
                ax.set_xlabel(u'Posição (µm)')
                ax.set_ylabel('Teor de carbono (% massa)')

        for ax in axes[i+1:]:
            ax.set_axis_off()

        ax = axes[0]
        cmax = x2wp(cmax, y=y)
        ax.annotate(s='{:.2f} %'.format(cmax), xy=(zmax, cmax),
                    xytext=(12, -10), textcoords='offset points', ha='left',
                    arrowprops=dict(arrowstyle='->'))

        if args.save:
            fname = os.path.join(args.dir, cprofiles.basename + '_sep.svg')
            fig.savefig(fname, bbox_inches='tight')
            os.system('svg2pdf ' + fname)

    if args.silent:
        plt.show()
    # plt.close()
