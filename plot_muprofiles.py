# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt
from cpartition import FCC, BCC, Interface, WBs, CProfiles, x2wp

if __name__ == '__main__':
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
    intf = Interface(domain1=aust, domain2=ferr, type_int='mobile.eq')
    # WBs composition
    _, cwbs = intf.comp()
    cwbs = x2wp(cwbs, y=y)

    basename = 'coupled_FoFo_375_CCE'
    mumart = None
    # tlist = [0.1, 1, 10, 60, 100, 1000]
    # basename = 'coupled_FoFo_375_mu23e3'
    # mumart = 23207
    # basename = 'coupled_FoFo_375_mu20e3'
    # mumart = 20e3
    # basename = 'coupled_FoFo_375_mu30e3'
    # mumart = 30e3
    # basename = 'coupled_FoFo_375_CCEpara'
    # mumart = 35511.1
    csolubility = 5e-4
    tlist = [0.1, 1, 10, 100, 1000]

    try:
        tracking
    except:
        cprofiles = CProfiles(basename, 'C_profiles')
        
    fig, ax = plt.subplots(figsize=(6, 4))

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

    for t in tlist:
        j, = cprofiles.where_tlist([t], [])
        
        strct = cprofiles.ss[j]
        cprofiles.plot_cprofiles(ax=ax, mirror=True,
                                 func=lambda x: x2mu(x, strct, mumart, csolubility),
                                 tlist=[t])

    ax.set_xlim(-1.16, 1.16)
    # ax.set_ylim(10, 60)
    ax.set_xlabel(u'Posição (µm)')
    ax.set_ylabel(r'$\mu_C$ (kJ/mol)')
    ax.legend(fancybox=False)

    fname = '/home/arthur/tese/img/cpartition/muprofiles/' + basename + '.svg'
    plt.savefig(fname, bbox_inches='tight')
    import os
    os.system('svg2pdf ' + fname)

    # plt.show()
    # plt.close()
