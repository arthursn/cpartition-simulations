#!/usr/bin/python3

# -*- coding: utf-8 -*-

if __name__ == '__main__':
    import os
    import numpy as np
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    from cpartition import FCC, BCC, WBs, CProfiles
    import argparse

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

        sel = (strct == 'fer1') | (strct == 'fer2') | (strct == 'fer3')
        mu[sel] = ferr.x2mu['C'](c[sel])

        sel = (strct == 'aus1') | (strct == 'aus2') | (strct == 'aust')
        mu[sel] = aust.x2mu['C'](c[sel])

        return 1e-3*mu

    dictloc = dict(ur='upper right', ul='upper left',
                   ll='lower left', lr='lower right',
                   best='best')

    parser = argparse.ArgumentParser()
    parser.add_argument('basenames', nargs='+')
    parser.add_argument('-s', '--silent', action='store_false')
    parser.add_argument('-S', '--save', action='store_true')

    parser.add_argument('-d', '--dir', default='')
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

    parser.add_argument('-M', '--mu', type=float, default=None)
    parser.add_argument('-c', '--sol', type=float, default=5.4e-4)

    args = parser.parse_args()

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

        try:
            fig, ax = plt.subplots(figsize=args.figsize)
            cprofiles = CProfiles(basename, 'C_profiles')

            last_t = None
            for t in args.time:
                j = cprofiles.where_tlist([t], [])
                if len(j) > 0:
                    j = j[0]
                    strct = cprofiles.ss[j]
                    cprofiles.plot_cprofiles(ax=ax, mirror=args.mirror,
                                             func=lambda x: x2mu(
                                                 x, strct, args.mu,
                                                 args.sol),
                                             tlist=[t])
                    last_t = t
        except Exception:
            print('Failed to plot "{}"'.format(basename))
            plt.close()
        else:
            ax.set_xlim(args.xlim)
            ax.set_ylim(args.ylim)
            ax.set_xlabel(u'Posição (µm)')
            ax.set_ylabel(r'$\mu_C$ (kJ/mol)')
            ax.legend(loc=dictloc[args.loc], fancybox=False)

            if args.tracking:
                cprofiles.plot_locus_interface([('aus1.s0', 'aus1.ci0'),
                                                ('aus1.sn', 'aus1.cin'),
                                                ('aus2.s0', 'aus2.ci0'),
                                                ('aus2.sn', 'aus2.cin')],
                                               ax=ax, mirror=args.mirror,
                                               func=lambda x: 1e-3 *
                                               aust.x2mu['C'](x),
                                               color='k', ls=':',
                                               lw=.8, label='')

            if args.label:
                if last_t is not None:
                    cprofiles.label_phases(ax=ax, t=last_t,
                                           labels=labels,
                                           mirror=args.mirror, size=12)

            if args.save:
                fname = os.path.join(args.dir,
                                     cprofiles.basename +
                                     args.append + args.ext)
                plt.savefig(fname, bbox_inches='tight')
                os.system('svg2pdf ' + fname)

    if args.silent:
        plt.show()
    else:
        plt.close('all')
