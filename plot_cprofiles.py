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
    from parse_args import lookup_option, split_string

    rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': 'Arial',
                     'font.size': 13})

    y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
             Si=5.02504411E-2, Fe=9.4414085022e-1)

    dictloc = dict(ur='upper right', ul='upper left',
                   ll='lower left', lr='lower right',
                   best='best')

    args = sys.argv[1:]

    if len(args) > 0:
        # Saving options
        save, args = lookup_option('-save', args, None, [])
        save = True if len(save) > 0 else False

        directory, args = lookup_option('-dir', args, str, [])
        directory = directory[-1] if len(
            directory) > 0 else '/home/arthur/tese/img/cpartition/cprofiles/'

        suffix, args = lookup_option('-suffix', args, str, [])
        suffix = suffix[-1] if len(suffix) > 0 else ''

        ext, args = lookup_option('-ext', args, str, [])
        ext = '.' + ext[-1].strip('.') if len(ext) > 0 else '.svg'

        # Plotting options
        title, args = lookup_option('-title', args, None, [])
        title = True if len(title) > 0 else False

        log, args = lookup_option('-log', args, None, [])
        log = True if len(log) > 0 else False

        loc, args = lookup_option('-loc', args, str, [])
        loc = dictloc[loc[-1]] if len(loc) > 0 else 'best'

        xlim, args = lookup_option('-xlim', args, str, [])
        xlim = split_string(xlim[-1], float) if len(xlim) > 0 else (None, None)

        ylim, args = lookup_option('-ylim', args, str, [])
        ylim = split_string(ylim[-1], float) if len(ylim) > 0 else (None, None)

        show, args = lookup_option('-show', args, None, [])
        show = True if len(show) > 0 else False

        figsize, args = lookup_option('-figsize', args, float, [], True)
        figsize = figsize if len(figsize) == 2 else (6, 4)

        # Options passed to cprofiles
        label, args = lookup_option('-label', args, None, [])
        label = True if len(label) > 0 else False

        mirror, args = lookup_option('-mirror', args, None, [])
        mirror = True if len(mirror) > 0 else False

        tracking, args = lookup_option('-tracking', args, None, [])
        tracking = True if len(tracking) > 0 else False

        tlist, args = lookup_option('-t', args, float, [], True)
        tlist = sorted(tlist)

        every, args = lookup_option('-every', args, str, [])
        tmin, tmax, each = split_string(
            every[-1]) if len(every) > 0 else (None, None, None)

        for basename in args:
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
                fig, ax = plt.subplots(figsize=figsize)
                cprofiles = CProfiles(basename)
                cprofiles.plot_cprofiles(ax=ax, slc=slice(tmin, tmax, each),
                                         tlist=tlist, mirror=mirror,
                                         func=lambda x: x2wp(x, y=y))

            except:
                print('Failed to plot "{}"'.format(basename))
                plt.close(fig)
            else:
                ax.set_xlim(*xlim)
                ax.set_ylim(*ylim)
                ax.set_xlabel(u'Posição (μm)')
                ax.set_ylabel('Teor de carbono (%)')
                ax.legend(loc=loc, fancybox=False)

                if tracking:
                    cprofiles.plot_locus_interface(tracked_interfaces, ax=ax, mirror=mirror,
                                                   color='k', ls='--', lw=.8,
                                                   func=lambda x: x2wp(x, y=y), label='')

                if log:
                    ax.set_yscale('log')

                if title:
                    ax.set_title(cprofiles.basename)

                if label:
                    cprofiles.label_phases(ax=ax, t=tlist[-1],
                                           labels=labels,
                                           mirror=mirror, size=12)

                if save:
                    fname = os.path.join(
                        directory, cprofiles.basename + suffix + ext)
                    plt.savefig(fname, bbox_inches='tight')
                    os.system('svg2pdf ' + fname)

        if show:
            plt.show()
        else:
            plt.close('all')
    else:
        print('Nothing to plot')
