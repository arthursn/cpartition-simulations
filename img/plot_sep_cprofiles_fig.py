# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cpartition import FCC, BCC, Interface, WBs, CProfiles, x2wp

y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
         Si=5.02504411E-2, Fe=9.4414085022e-1)

# For calculating WBs
T_C = 375.
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
    pos = pd.read_table('../pos_extremities/{}.txt'.format(basename), sep=' ')
    ci = pd.read_table('../C_extremities/{}.txt'.format(basename), sep=' ')

tracking = False
ncol = 2
nrow = len(tlist)//ncol

fig, axes = plt.subplots(nrow, ncol, figsize=(4*ncol, 3*nrow))
axes = axes.ravel()

for i, (t, ax) in enumerate(zip(tlist, axes)):
    cprofiles.plot_cprofiles(ax=ax, mirror=True,
                             func=lambda x: x2wp(x, y=y),
                             tlist=[t], color='k', lw=1)
    if i == 0:
        j, = cprofiles.where_tlist([t], [])
        idxmax = cprofiles.cc[j].argmax()
        zmax = cprofiles.zz[j][idxmax]
        cmax = cprofiles.cc[j][idxmax]

    if tracking:
        ax.plot(pos['aus1.sn'], x2wp(ci['aus1.cin'], y=y), 'k:')

    ax.axhline(cwbs, color='k', ls='-.')
    ax.text(-1.1, cwbs + .1, s='WBs', ha='left')
    ax.text(.05, .9, s='{} s'.format(t), transform=ax.transAxes)

    ax.set_xlim(-1.16, 1.16)
    if i > 2:
        ax.set_ylim(0, 2.5)
    else:
        ax.set_ylim(0)

    if i == len(tlist) - ncol:
        ax.set_xlabel(u'Posição (µm)')
        ax.set_ylabel('Teor de carbono (% massa)')

ax = axes[0]
cmax = x2wp(cmax, y=y)
ax.annotate(s='{:.2f} %'.format(cmax), xy=(zmax, cmax),
            xytext=(20, -10), textcoords='offset points',
            arrowprops=dict(arrowstyle='->'))

plt.savefig(('/home/arthur/tese/img/cpartition/'
             'cprofiles_sep_coupled_CCE_.svg'), bbox_inches='tight')
plt.show()
# plt.close()
