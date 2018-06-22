# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt

# script in the local folder
from plot_cavg import plot_cavg

from matplotlib import rcParams
rcParams.update({'font.family': 'sans-serif',
                 'font.sans-serif': 'Arial',
                 'font.size': 13})

files = ['../C_avg/coupled_FoFo_375_CCEortho.txt',
         '../C_avg/coupled_FoFo_375_mu20e3.txt',
         '../C_avg/coupled_FoFo_375_mu23e3.txt',
         '../C_avg/coupled_FoFo_375_mu30e3.txt',
         '../C_avg/coupled_FoFo_375_CCEpara.txt']
labels = [r"$\alpha'-\theta$" + u' ortho',
          r'$\mu_C =$' + u' 20.0 kJ/mol',
          r'$\mu_C =$' + u' 23.2 kJ/mol',
          r'$\mu_C =$' + u' 30.0 kJ/mol',
          r"$\alpha'-\theta$" + u' para']

fig, ax = plt.subplots(figsize=(6, 4))

plot_cavg(files, labels, ax=ax, lw=1)
ax.axhline(c0, ls=':', color='k', lw=1)
ax.axhline(WBs, ls=':', color='k', lw=1)
add_label(ax, 'd)', py=0)

fout = 'cpartition.pdf'
fig.savefig(fout, bbox_inches='tight')

plt.show()
plt.close('all')
