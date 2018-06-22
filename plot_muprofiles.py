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

        sel = (strct == 'aus1') | (strct == 'aus2')
        mu[sel] = aust.x2mu['C'](c[sel])

        return 1e-3*mu

    # basename = 'coupled_FoFo_375_CCE'
    # mumart = None
    # tlist = [0.1, 1, 10, 60, 100, 1000]
    # basename = 'coupled_FoFo_375_mu23e3'
    # mumart = 23207
    # basename = 'coupled_FoFo_375_mu20e3'
    # mumart = 20e3
    # basename = 'coupled_FoFo_375_mu30e3'
    # mumart = 30e3
    # basename = 'coupled_FoFo_375_CCEortho'
    # mumart = 8188.68
    # basename = 'coupled_FoFo_375_CCEpara'
    # mumart = 35511.1
    # csolubility = 5e-4
    # tlist = [0.1, 1, 10, 100, 1000]

    dictloc = dict(ur='upper right', ul='upper left',
                   ll='lower left', lr='lower right',
                   best='best')

    args = sys.argv[1:]

    if len(args) > 0:
        # Saving options
        save, args = lookup_option('-save', args, None, [])
        save = True if len(save) > 0 else False

        directory, args = lookup_option('-dir', args, str, [])
        directory = directory[-1] if len(directory) > 0 else '/home/arthur/tese/img/cpartition/muprofiles/'

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
        tlist, args = lookup_option('-t', args, float, [], True)
        tlist = sorted(tlist)

        mumart, args = lookup_option('-mu', args, float, [])
        mumart = mumart[-1] if len(mumart) > 0 else None 
        
        csolubility, args = lookup_option('-sol', args, float, [])
        csolubility = csolubility[-1] if len(csolubility) > 0 else 5.4e-4

        for basename in args:
            fig, ax = plt.subplots(figsize=figsize)

            cprofiles = CProfiles(basename, 'C_profiles')
            
            for t in tlist:
                j, = cprofiles.where_tlist([t], [])
                
                strct = cprofiles.ss[j]
                cprofiles.plot_cprofiles(ax=ax, mirror=True,
                                         func=lambda x: x2mu(x, strct, mumart, csolubility),
                                         tlist=[t])

            ax.set_xlim(*xlim)
            ax.set_xlim(*ylim)

            ax.set_xlabel(u'Posição (µm)')
            ax.set_ylabel(r'$\mu_C$ (kJ/mol)')
            ax.legend(fancybox=False)

            if save:
                fname = os.path.join(directory, cprofiles.basename + suffix + ext)
                plt.savefig(fname, bbox_inches='tight')
                os.system('svg2pdf ' + fname)

        if show:
            plt.show()
        else:
            plt.close('all')
    else:
        print('Nothing to plot')
