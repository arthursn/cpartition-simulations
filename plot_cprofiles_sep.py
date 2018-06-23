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
    from parse_args import lookup_option, split_string

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

    args = sys.argv[1:]

    if len(args) > 1:
        # Saving options
        save, args = lookup_option('-save', args, None, [])
        save = True if len(save) > 0 else False

        directory, args = lookup_option('-dir', args, str, [])
        directory = directory[-1] if len(directory) > 0 else '/home/arthur/tese/img/cpartition/cprofiles/'

        suffix, args = lookup_option('-suffix', args, str, [])
        suffix = suffix[-1] if len(suffix) > 0 else ''

        ext, args = lookup_option('-ext', args, str, [])
        ext = '.' + ext[-1].strip('.') if len(ext) > 0 else '.svg'

        # Plotting options
        xlim, args = lookup_option('-xlim', args, str, [])
        xlim = split_string(xlim[-1], float) if len(xlim) > 0 else (None, None)

        ylim, args = lookup_option('-ylim', args, str, [])
        ylim = split_string(ylim[-1], float) if len(ylim) > 0 else (None, None)

        show, args = lookup_option('-show', args, None, [])
        show = True if len(show) > 0 else False

        figsize, args = lookup_option('-figsize', args, float, [], True)
        figsize = figsize if len(figsize) == 2 else (6, 4)

        # Options passed to cprofiles
        mirror, args = lookup_option('-mirror', args, None, [])
        mirror = True if len(mirror) > 0 else False

        tracking, args = lookup_option('-tracking', args, None, [])
        tracking = True if len(tracking) > 0 else False

        tlist, args = lookup_option('-t', args, float, [], True)
        tlist = sorted(tlist)

        # Special bois
        mommysaysimspecial, args = lookup_option('-special', args, None, [])
        mommysaysimspecial = True if len(mommysaysimspecial) > 0 else False

        wbs, args = lookup_option('-wbs', args, None, [])
        wbs = True if len(wbs) > 0 else False

        miniall, args = lookup_option('-all', args, None, [])
        miniall = True if len(miniall) > 0 else False

        if miniall:
            tlist.append(tlist.copy())

        ncol = 2
        nrow = int(np.ceil(len(tlist)/ncol))

        for basename in args:
            cprofiles = CProfiles(basename, 'C_profiles')

            if tracking:
                pos = pd.read_table(
                    'pos_extremities/{}.txt'.format(cprofiles.basename), sep=' ')
                ci = pd.read_table(
                    'C_extremities/{}.txt'.format(cprofiles.basename), sep=' ')

            fig, axes = plt.subplots(nrow, ncol, figsize=(4*ncol, 3*nrow))
            # plt.subplots_adjust(wspace=.5, hspace=.5)
            axes = axes.ravel()

            for i, (t, ax) in enumerate(zip(tlist, axes)):
                kw = dict(lw=1)
                if not isinstance(t, list):
                    t = [t]
                    kw['color'] = next(colorcycle)

                cprofiles.plot_cprofiles(ax=ax, mirror=True,
                                         func=lambda x: x2wp(x, y=y),
                                         tlist=t, **kw)
                if i == 0:
                    j, = cprofiles.where_tlist(t, [])
                    idx, = np.where(cprofiles.ss[j] == 'aus1')
                    idx = idx[0]
                    zmax = 2*cprofiles.zz[j][-1] - cprofiles.zz[j][idx]
                    cmax = cprofiles.cc[j][idx]

                if tracking:
                    ax.plot(pos['aus1.sn'], x2wp(ci['aus1.cin'], y=y), 'k:')

                ax.set_ylim(*ylim)
                ax.set_xlim(*xlim)

                if len(t) == 1:
                    ax.text(.01, .85, s='{:g} s'.format(t[0]), transform=ax.transAxes)

                    
                    if wbs:
                        ax.axhline(cwbs, color='k', ls='--', lw=.8)
                        ax.text(1.1, cwbs + .1, s='WBs', ha='right')

                    # special axis dimensions
                    if mommysaysimspecial:
                        if i > 1:
                            ax.set_ylim(-.1, 2.5)
                        else:
                            ax.set_ylim(-.1)
                    
                    cprofiles.label_phases(ax, t, mirror=True)
                else:
                    ax.text(.01, .9, s='Todos', transform=ax.transAxes)

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

            if save:
                fname = os.path.join(
                    directory, cprofiles.basename + '_sep.svg')
                fig.savefig(fname, bbox_inches='tight')
                os.system('svg2pdf ' + fname)

        if show:
            plt.show()
        # plt.close()
