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

files = ['../C_avg/mart_FoFo_CCEortho.txt',
         '../C_avg/mart_FoFo_mu16e3.txt',
         '../C_avg/mart_FoFo_mu20e3.txt',
         '../C_avg/mart_FoFo_mu23e3.txt',
         '../C_avg/mart_FoFo_mu25e3.txt',
         '../C_avg/mart_FoFo_mu30e3.txt',
         '../C_avg/mart_FoFo_CCEpara.txt',
         '../C_avg/mart_FoFo_CCE.txt']
labels = [r"$\alpha'-\theta$" + u' ortho',
          r'$\mu_C =$' + u' 16,4 kJ/mol',
          r'$\mu_C =$' + u' 20,0 kJ/mol',
          r'$\mu_C =$' + u' 23,2 kJ/mol',
          r'$\mu_C =$' + u' 25,0 kJ/mol',
          r'$\mu_C =$' + u' 30,0 kJ/mol',
          r"$\alpha'-\theta$" + u' para',
          'Sem carbonetos']

fig, ax = plt.subplots(figsize=(7, 5))

y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
         Si=5.02504411E-2, Fe=9.4414085022e-1)
c0 = x2wp(3.34414e-02, y=y)
WBs = 1.5902

ax = plot_cavg(files, labels, ax=ax, lw=1, y=y, plotcarbides=False)

ax.axhline(c0, ls=':', color='k', lw=1)
ax.axhline(WBs, ls=':', color='k', lw=1)

ax.set_xscale('log')
ax.set_xlim(1e-2, 1e4)
ax.set_ylim(0, 1.4)
ax.set_xlabel('Tempo (s)')
ax.set_ylabel(u'Carbono médio em ' + r"$\alpha' + \theta$ (%)")
# ax2.set_ylabel(u'Fração estimada de ' + r'$\theta$ (% vol.)' )
ax.legend(loc='best', fancybox=False, ncol=2)

fout = '/home/arthur/tese/img/cpartition/nobainite_cavg.svg'
fig.savefig(fout, bbox_inches='tight')
os.system('svg2pdf ' + fout)

plt.show()
# plt.close('all')
