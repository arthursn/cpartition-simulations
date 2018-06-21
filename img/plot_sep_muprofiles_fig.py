# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cpartition import FCC, BCC, Interface, WBs, CProfiles, x2wp

y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
         Si=5.02504411E-2, Fe=9.4414085022e-1)

# For calculating WBs
T_C = 375.
mart = BCC(T_C=T_C, tdata='../thermo/FoFo/TCFE8/375-bcc.txt')
aust = FCC(T_C=T_C, tdata='../thermo/FoFo/TCFE8/375-fcc.txt')
ferr = BCC(T_C=T_C, tdata='../thermo/FoFo/TCFE8/375-bcc.txt', E=WBs(T_C))
intf = Interface(domain1=aust, domain2=ferr, type_int='mobile.eq')
# WBs composition
_, cwbs = intf.comp()
cwbs = x2wp(cwbs, y=y)

basename = 'coupled_FoFo_375_CCE'
tlist = [0.1, 1, 10, 60, 100, 1000]

try:
    tracking
except:
    cprofiles = CProfiles(basename, '../C_profiles')

tracking = False
ncol = 2
nrow = len(tlist)//ncol

fig, axes = plt.subplots(nrow, ncol, figsize=(4*ncol, 3*nrow))
axes = axes.ravel()

for i, (t, ax) in enumerate(zip(tlist, axes)):
    j, = cprofiles.where_tlist([t], [])
    
    strct = cprofiles.ss[j]
    def x2mu(c):
        mu = np.full(c.shape, 0)
        
        sel = (strct == 'mart')
        mu[sel] = mart.x2mu['C'](c[sel])
        
        sel = (strct == 'fer1') | (strct == 'fer2')
        mu[sel] = ferr.x2mu['C'](c[sel])

        sel = (strct == 'aus1') | (strct == 'aus2')
        mu[sel] = aust.x2mu['C'](c[sel])

        return 1e-3*mu

    cprofiles.plot_cprofiles(ax=ax, mirror=True,
                             func=lambda x: x2mu(x),
                             tlist=[t], color='k', lw=1)

    ax.set_xlim(-1.16, 1.16)
    ax.set_ylim(10, 60)

    if i == len(tlist) - ncol:
        ax.set_xlabel(u'Posição (µm)')
        ax.set_ylabel(r'$\mu_C$ (kJ/mol)')

# plt.savefig(('/home/arthur/tese/img/cpartition/'
             # 'muprofiles_coupled_CCE_.svg'), bbox_inches='tight')
plt.show()
# plt.close()
