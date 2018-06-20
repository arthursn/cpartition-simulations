#!/usr/bin/python3
# -*- coding: utf-8 -*-


def lookup_option(option, args, dtype=None, vallist=[]):
    if option in args:
        idx = args.index(option)
        args.pop(idx)
        if dtype:
            try:
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


def parse_slice(slc, dtype=int):
    spt = []
    try:
        spt = list(map(lambda s: filter_number(s, dtype), slc.split(':')))
    except:
        print('Failed parsing {}'.format(slc))
    return spt


if __name__ == '__main__':
    import sys
    import os
    import numpy as np
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    from cpartition import x2wp, CProfiles

    rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': 'Arial',
                     'font.size': 13,
                     'mathtext.fontset': 'stix'})

    y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
             Si=5.02504411E-2, Fe=9.4414085022e-1)

    args = sys.argv[1:]

    if len(args) > 0:
        mirror = False
        mirror, args = lookup_option('-m', args, None, [])
        if len(mirror) > 0:
            mirror = True
        
        log = False
        log, args = lookup_option('-l', args, None, [])
        if len(log) > 0:
            log = True

        sel = []
        sel, args = lookup_option('-s', args, int, [])
        sel = sorted(sel)

        t_sel = []
        t_sel, args = lookup_option('-t', args, float, [])
        t_sel = sorted(t_sel)

        tmin, tmax, each = None, None, None
        every, args = lookup_option('-e', args, str, [])
        if len(every) > 0:
            tmin, tmax, each = parse_slice(every[-1])

        xlim = (None, None)
        xlim, args = lookup_option('-xlim', args, str, [])
        if len(xlim) > 0:
            xlim = parse_slice(xlim[-1], float)

        ylim = (None, None)
        ylim, args = lookup_option('-ylim', args, str, [])
        if len(ylim) > 0:
            ylim = parse_slice(ylim[-1], float)

        for basename in args:
            try:
                fig, ax = plt.subplots(figsize=(8, 5))
                cprofiles = CProfiles(basename)
                cprofiles.load_time()
                
                if len(t_sel) > 0:
                    for t in t_sel:
                        matches, = np.where(cprofiles.t == t)
                        sel += list(matches)

                cprofiles.plot_cprofiles(ax=ax,
                                        slc=slice(tmin, tmax, each), sel=sel,
                                        mirror=mirror, func=lambda x: x2wp(x, y=y),
                                        vmin=0, vmax=1.8)
            except:
                raise
                print('Failed to plot "{}"'.format(basename))
                plt.close(fig)
            else:
                ax.set_xlabel(u'Position (μm)')
                ax.set_ylabel('Carbon content (wt.%)')
                ax.legend(loc='upper left', fancybox=False)
                ax.set_xlim(*xlim)
                ax.set_ylim(*ylim)
                if log:
                    ax.set_yscale('log')
                ax.set_title(basename)

            # plt.savefig('img/' + basename + '.png', dpi=150)
            # plt.close()

        plt.show()
    else:
        print('Nothing to plot')
