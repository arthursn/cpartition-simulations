#!/usr/bin/python3
# -*- coding: utf-8 -*-


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
        spt = list(map(lambda s: filter_number(s, dtype), string.split(splitchar)))
    except:
        print('Failed parsing {}'.format(string))
    return spt


if __name__ == '__main__':
    import sys
    import os
    import numpy as np
    from fnmatch import fnmatch
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    from cpartition import x2wp, CProfiles

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

        title, args = lookup_option('-title', args, None, [])
        if len(title) > 0:
            title = True
        else:
            title = False

        directory, args = lookup_option('-dir', args, str, [])
        if len(directory) > 0:
            directory = directory[-1]
        else:
            directory = '/home/arthur/tese/img/cpartition/cprofiles/'

        mirror, args = lookup_option('-mirror', args, None, [])
        if len(mirror) > 0:
            mirror = True
        else:
            mirror = False

        log, args = lookup_option('-log', args, None, [])
        if len(log) > 0:
            log = True
        else:
            log = False

        loc, args = lookup_option('-loc', args, str, [])
        if len(loc) > 0:
            loc = dictloc[loc[-1]]
        else:
            loc = 'best'

        tracking, args = lookup_option('-tracking', args, None, [])
        if len(tracking) > 0:
            tracking = True
        else:
            tracking = False

        tlist, args = lookup_option('-t', args, float, [], True)
        tlist = sorted(tlist)

        figsize, args = lookup_option('-figsize', args, float, [], True)
        if len(figsize) != 2:
            figsize = [6, 4]

        every, args = lookup_option('-every', args, str, [])
        if len(every) > 0:
            tmin, tmax, each = split_string(every[-1])
        else:
            tmin, tmax, each = None, None, None

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

        suffix, args = lookup_option('-suffix', args, str, [])
        if len(suffix) > 0:
            suffix = suffix[0]
        else:
            suffix = ''

        ext, args = lookup_option('-ext', args, str, [])
        if len(ext) > 0:
            ext = '.' + ext[0].strip('.')
        else:
            ext = '.svg'

        for basename in args:
            try:
                fig, ax = plt.subplots(figsize=figsize)
                cprofiles = CProfiles(basename)
                cprofiles.plot_cprofiles(ax=ax, slc=slice(tmin, tmax, each),
                                         tlist=tlist, mirror=mirror,
                                         func=lambda x: x2wp(x, y=y))

                if tracking:
                    cprofiles.plot_locus_interface([('aus1.sn', 'aus1.cin'),
                                                    ('aus2.s0', 'aus2.ci0'),
                                                    ('aus2.sn', 'aus2.cin')],
                                                   ax=ax, color='k', ls='--', lw=1,
                                                   func=lambda x: x2wp(x, y=y), label='')
            except:
                raise
                print('Failed to plot "{}"'.format(basename))
                plt.close(fig)
            else:
                ax.set_xlabel(u'Posição (μm)')
                ax.set_ylabel('Teor de carbono (%)')
                ax.legend(loc=loc, fancybox=False)
                ax.set_xlim(*xlim)
                ax.set_ylim(*ylim)
                if log:
                    ax.set_yscale('log')

                if title:
                    ax.set_title(cprofiles.basename)

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
