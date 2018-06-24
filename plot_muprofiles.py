# -*- coding: utf-8 -*-

if __name__ == '__main__':
    import sys
    import os
    import numpy as np
    import pandas as pd
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    from cpartition import FCC, BCC, Interface, WBs, CProfiles, x2wp
    from parse_args import lookup_option, split_string

    rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': 'Arial',
                     'font.size': 13})

    y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
             Si=5.02504411E-2, Fe=9.4414085022e-1)

    # For calculating WBs
    T_C = 375.
    mart = BCC(T_C=T_C, tdata='thermo/FoFo/TCFE8/375-bcc.txt')
    aust = FCC(T_C=T_C, tdata='thermo/FoFo/TCFE8/375-fcc.txt')
    ferr = BCC(T_C=T_C, tdata='thermo/FoFo/TCFE8/375-bcc.txt', E=WBs(T_C))

    def x2mu(c, strct, mumart=None, csolubility=5.4e-4):
        mu = np.full(c.shape, 0)
        sel = (strct == 'mart')
        if mumart:
            # CCEtheta case
            selbyc = sel & (c > csolubility)
            mu[selbyc] = mumart
            selbyc = sel & (c <= csolubility)
            mu[selbyc] = mart.x2mu['C'](c[selbyc])
        else:
            mu[sel] = mart.x2mu['C'](c[sel])

        sel = (strct == 'fer1') | (strct == 'fer2')
        mu[sel] = ferr.x2mu['C'](c[sel])

        sel = (strct == 'aus1') | (strct == 'aus2') | (strct == 'aust')
        mu[sel] = aust.x2mu['C'](c[sel])

        return 1e-3*mu

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
            directory) > 0 else '/home/arthur/tese/img/cpartition/muprofiles/'

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
        label, args = lookup_option('-label', args, None, [])
        label = True if len(label) > 0 else False

        tracking, args = lookup_option('-tracking', args, None, [])
        tracking = True if len(tracking) > 0 else False

        tlist, args = lookup_option('-t', args, float, [], True)
        tlist = sorted(tlist)

        mumart, args = lookup_option('-mu', args, float, [])
        mumart = mumart[-1] if len(mumart) > 0 else None

        csolubility, args = lookup_option('-sol', args, float, [])
        csolubility = csolubility[-1] if len(csolubility) > 0 else 5.4e-4

        for basename in args:
            if 'mart' in basename:
                labels = [('aus1', r'$\gamma_1$', 1),
                          ('aus2', r'$\gamma_2$', 1),
                          ('aust', r'$\gamma$', 1),
                          ('mart', r"$\alpha'$", 1)]
            else:
                labels = [('aus1', r'$\gamma_1$', -1),
                          ('aus2', r'$\gamma_2$', -1),
                          ('aust', r'$\gamma$', -1),
                          ('mart', r"$\alpha'$", -1),
                          ('fer1', r'$\alpha_{b1}$', 1),
                          ('fer2', r'$\alpha_{b2}$', 1)]

            try:
                fig, ax = plt.subplots(figsize=figsize)
                cprofiles = CProfiles(basename, 'C_profiles')

                for t in tlist:
                    j, = cprofiles.where_tlist([t], [])

                    strct = cprofiles.ss[j]
                    cprofiles.plot_cprofiles(ax=ax, mirror=True,
                                             func=lambda x: x2mu(x, strct, mumart, csolubility),
                                             tlist=[t])
            except:
                print('Failed to plot "{}"'.format(basename))
                plt.close(fig)
            else:
                ax.set_xlim(*xlim)
                ax.set_ylim(*ylim)
                ax.set_xlabel(u'Posição (µm)')
                ax.set_ylabel(r'$\mu_C$ (kJ/mol)')
                ax.legend(fancybox=False)

                if tracking:
                    cprofiles.plot_locus_interface([('aus1.s0', 'aus1.ci0'),
                                                    ('aus1.sn', 'aus1.cin'),
                                                    ('aus2.s0', 'aus2.ci0'),
                                                    ('aus2.sn', 'aus2.cin')],
                                                   ax=ax, mirror=True,
                                                   func=lambda x: 1e-3*aust.x2mu['C'](x),
                                                   color='k', ls=':', lw=.8, label='')

                if label:
                    cprofiles.label_phases(ax=ax, t=tlist[-1],
                                           labels=labels,
                                           mirror=True, size=12)


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
