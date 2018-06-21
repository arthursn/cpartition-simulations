# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cpartition import FCC, BCC, Interface, WBs, CProfiles, x2wp


def lookup_option(option, args, dtype=None, vallist=[], multi=False):
    """
    Get options in list of arguments
    """
    if option in args:
        idx = args.index(option)
        args.pop(idx)
        if dtype:
            try:
                if multi:
                    while True:
                        try:
                            val = args[idx]
                            val = dtype(val)
                        except:
                            break
                        else:
                            vallist += [val]
                            args.pop(idx)
                else:
                    val = args.pop(idx)
                    val = dtype(val)
                    vallist += [val]
            except IndexError:
                print('No argument provided')
            except ValueError:
                print('Failed parsing {} as {}'.format(val, dtype))
            except:
                print('Unexpected error')
        else:
            vallist += [True]

    while option in args:
        vallist, args = lookup_option(option, args, dtype, vallist)

    return vallist, args


def filter_number(s, dtype=int):
    if s == '':
        s = None
    else:
        try:
            s = dtype(s)
        except:
            raise
    return s


def split_string(string, dtype=int, splitchar=':'):
    spt = []
    try:
        spt = list(map(lambda s: filter_number(
            s, dtype), string.split(splitchar)))
    except:
        print('Failed parsing {}'.format(string))
    return spt


if __name__ == '__main__':
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
        show, args = lookup_option('-show', args, None, [])
        if len(show) > 0:
            show = True
        else:
            show = False

        save, args = lookup_option('-save', args, None, [])
        if len(save) > 0:
            save = True
        else:
            save = False

        directory, args = lookup_option('-dir', args, str, [])
        if len(directory) > 0:
            directory = directory[-1]
        else:
            directory = '/home/arthur/tese/img/cpartition/cprofiles/'

        tracking, args = lookup_option('-tracking', args, None, [])
        if len(tracking) > 0:
            tracking = True
        else:
            tracking = False

        tlist, args = lookup_option('-t', args, float, [], True)
        tlist = sorted(tlist)

        xlim, args = lookup_option('-xlim', args, str, [])
        if len(xlim) > 0:
            xlim = split_string(xlim[-1], float)
        else:
            xlim = (None, None)

        ylim, args = lookup_option('-ylim', args, str, [])
        if len(ylim) > 0:
            ylim = split_string(ylim[-1], float)
        else:
            ylim = (None, None)

        ncol = 2
        nrow = len(tlist)//ncol
        
        for basename in args:
            cprofiles = CProfiles(basename, 'C_profiles')
            
            if tracking:
                pos = pd.read_table(
                    'pos_extremities/{}.txt'.format(cprofiles.basename), sep=' ')
                ci = pd.read_table(
                    'C_extremities/{}.txt'.format(cprofiles.basename), sep=' ')

            fig, axes = plt.subplots(nrow, ncol, figsize=(4*ncol, 3*nrow))
            axes = axes.ravel()

            for i, (t, ax) in enumerate(zip(tlist, axes)):
                cprofiles.plot_cprofiles(ax=ax, mirror=True,
                                         func=lambda x: x2wp(x, y=y),
                                         tlist=[t], color='k', lw=1)
                if i == 0:
                    j, = cprofiles.where_tlist([t], [])
                    idx, = np.where(cprofiles.ss[j] == 'aus1')
                    idx = idx[0]
                    zmax = 2*cprofiles.zz[j][-1] - cprofiles.zz[j][idx]
                    cmax = cprofiles.cc[j][idx]

                if tracking:
                    ax.plot(pos['aus1.sn'], x2wp(ci['aus1.cin'], y=y), 'k:')

                # ax.axhline(cwbs, color='k', ls='-.')
                # ax.text(1.1, cwbs + .1, s='WBs', ha='right')
                ax.text(.01, .9, s='{} s'.format(t), transform=ax.transAxes)

                ax.set_xlim(*xlim)
                ax.set_ylim(*ylim)
                # if i > 1:
                #     ax.set_ylim(0, 2.5)
                # else:
                #     ax.set_ylim(0)

                if i == len(tlist) - ncol:
                    ax.set_xlabel(u'Posição (µm)')
                    ax.set_ylabel('Teor de carbono (% massa)')

            ax = axes[0]
            cmax = x2wp(cmax, y=y)
            ax.annotate(s='{:.2f} %'.format(cmax), xy=(zmax, cmax),
                        xytext=(12, -10), textcoords='offset points', ha='left',
                        arrowprops=dict(arrowstyle='->'))

            if save:
                fname = os.path.join(directory, cprofiles.basename + '_sep.svg')
                fig.savefig(fname, bbox_inches='tight')
                os.system('svg2pdf ' + fname)

        if show:
            plt.show()
        # plt.close()
