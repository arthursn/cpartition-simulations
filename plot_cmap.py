#!/usr/bin/python3
# -*- coding: utf-8 -*-

if __name__ == '__main__':
    import os
    import sys
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    from cpartition import x2wp, CProfiles
    from parse_args import lookup_option

    rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': 'Arial',
                     'font.size': 13,
                     'mathtext.fontset': 'stix'})

    y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
             Si=5.02504411E-2, Fe=9.4414085022e-1)
    
    args = sys.argv[1:]
    
    if len(args) > 0:
        save, args = lookup_option('-save', args, None, [])
        save = True if len(save) > 0 else False

        for basename in args:
            cprofiles = CProfiles(basename)
            
            print(cprofiles.basename)

            ax = cprofiles.plot_colormap(mirror=True,
                                         func=lambda x: x2wp(x, y=y),
                                         vmin=0, vmax=1.8)
            ax.set_xlabel(u'Position (μm)')
            ax.set_ylabel('Time (s)')
            ax.set_title(cprofiles.basename)

            if save:
                plt.savefig(os.path.join('img', cprofiles.basename + '.png'), dpi=150)
                plt.close()

        plt.show()
    else:
        print('Nothing to plot')
