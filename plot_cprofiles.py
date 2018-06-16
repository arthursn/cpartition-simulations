#!/usr/bin/python3
# -*- coding: utf-8 -*-

def lookup_option(option, argv, dtype=None):
    val = None
    if option in argv:
        idx = argv.index(option)
        argv.pop(idx)
        try:
            val = argv.pop(idx)
            if dtype:
                val = dtype(val)
        except:
            print('Invalid option')

    return val, argv

if __name__ == '__main__':
    import sys
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    from cpartition import x2wp, CProfiles

    rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': 'Arial',
                     'font.size': 13,
                     'mathtext.fontset': 'stix'})

    y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
             Si=5.02504411E-2, Fe=9.4414085022e-1)

    if len(sys.argv) > 1:
        each, argv = lookup_option('-e', sys.argv, int)
        each, argv = lookup_option('--each', sys.argv, int)
        if not each:
            each = 100
        

        for fname in sys.argv[1:]:
            print(fname)

            cprofiles = CProfiles(fname)
            
            ax = cprofiles.plot_profiles(each=each, mirror=True,
                                         func=lambda x: x2wp(x, y=y),
                                         vmin=0, vmax=1.8)
            ax.set_xlabel(u'Position (μm)')
            ax.set_ylabel('Carbon content (wt.%)')
            ax.legend()

            ax.set_title(fname)
            
            # plt.savefig('img/' + fname.split('.')[0] + '.png', dpi=150)
            # plt.close()

        plt.show()
    else:
        print('Nothing to plot')
