# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt

from cpartition import x2wp

from matplotlib import rcParams
rcParams.update({'font.family': 'sans-serif',
                 'font.sans-serif': 'Arial',
                 'font.size': 13})

fname = '../C_avg/mart_FoFo_CCE.txt'
df = pd.read_csv(fname, sep=' ')

y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
         Si=5.02504411E-2, Fe=9.4414085022e-1)
c0 = x2wp(3.34414e-02, y=y)


fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(df['t'], x2wp(df['aust.cavg'], y=y), 'k--', label=r"$\gamma$")
ax.plot(df['t'], x2wp(df['mart.cavg'], y=y), 'k-', label=r"$\alpha'$")

ax.axhline(c0, ls=':', color='k', lw=1)
ax.set_xscale('log')
ax.set_xlim(1e-4, 1000)
# ax.set_ylim(0, 1.4)
ax.set_xlabel('Tempo (s)')
ax.set_ylabel(u'Teor de carbono médio na fase (% massa)')
ax.legend(loc='best', fancybox=False)
ax.annotate(s=r'$c_0$', xy=(.8e3, c0*1.05), ha='right')

fout = '/home/arthur/tese/img/cpartition/CCE_cavg.svg'
fig.savefig(fout, bbox_inches='tight')
os.system('svg2pdf ' + fout)

plt.show()
# plt.close('all')
