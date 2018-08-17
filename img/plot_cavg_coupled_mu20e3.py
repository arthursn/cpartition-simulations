# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt

from cpartition import x2wp

# script in the local folder
from plot_cavg import plot_cavg

from matplotlib import rcParams
rcParams.update({'font.family': 'sans-serif',
                 'font.sans-serif': 'Arial',
                 'font.size': 13})

files = ['../C_avg/coupled_FoFo_375_mu20e3.txt']
labels = [r'$\mu_C =$' + u' 20,0 kJ/mol']

fig, ax = plt.subplots(figsize=(5, 4))

y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
         Si=5.02504411E-2, Fe=9.4414085022e-1)
c0 = x2wp(3.34414e-02, y=y)
WBs = 1.5902

ax = plot_cavg(files, labels, ax=ax, lw=1, y=y, plotcarbides=False)

ax.axhline(c0, ls=':', color='k', lw=1)
ax.axhline(WBs, ls=':', color='k', lw=1)

ax.set_xscale('log')
ax.set_xlim(1e-2, 1000)
ax.set_ylim(0, 1.4)
ax.set_xlabel('Tempo (s)')
ax.set_ylabel(u'Carbono médio em ' + r"$\alpha' + \theta$ (%)")
ax.legend(loc='best', fancybox=False, ncol=2)


fout = '/home/arthur/tese/img/cpartition/coupled_cavg_mu20e3.svg'
fig.savefig(fout, bbox_inches='tight')
# os.system('svg2pdf ' + fout)

plt.show()
# plt.close('all')
