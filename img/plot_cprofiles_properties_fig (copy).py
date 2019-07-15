import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt

import os
from cpartition import x2wp, label, CProfiles

# script in the local folder
from plot_cavg import plot_cavg

rcParams.update({'font.family': 'sans-serif',
                 'font.sans-serif': 'Arial',
                 'font.size': 13,
                 'mathtext.fontset': 'stix'})

def add_label(ax, label, px=.15, py=.1, size=20, **kwargs):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x = float(min(xlim) - np.diff(xlim)*px)
    y = float(max(ylim) + np.diff(ylim)*py)
    ax.text(x=x, y=y, s=label, size=size, **kwargs)


def calculate_phase_fractions(fname):
    df = pd.read_csv(fname, delim_whitespace=True)

    domain_phase_map = dict()

    for colname in df.columns[1:]:
        domain = colname.split('.')[0]
        phase = ''.join([i for i in domain if not i.isdigit()])

        if phase not in domain_phase_map.keys():
            domain_phase_map[phase] = set([domain])
        else:
            domain_phase_map[phase].add(domain)

    total_length = df.iloc[0].max() - df.iloc[0].min()

    phase_fractions = pd.DataFrame(0, index=df.index, columns=domain_phase_map.keys())

    for phase, domains in domain_phase_map.items():
        for domain in domains:
            phase_fractions[phase] += (df['{}.sn'.format(domain)] - df['{}.s0'.format(domain)]).abs()

    phase_fractions = phase_fractions/total_length

    phase_fractions['t'] = df['t']
    cols = phase_fractions.columns.tolist()
    phase_fractions = phase_fractions[cols[-1:] + cols[:-1]]

    return phase_fractions


fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(9, 3.8))
plt.subplots_adjust(wspace=.3)

##########################################################

files = ['../pos_extremities/coupled_FoFo_375_CCEortho.txt',
         # '../pos_extremities/coupled_FoFo_375_mu16e3.txt',
         '../pos_extremities/coupled_FoFo_375_mu20e3.txt',
         '../pos_extremities/coupled_FoFo_375_mu23e3.txt',
         # '../pos_extremities/coupled_FoFo_375_mu25e3.txt',
         '../pos_extremities/coupled_FoFo_375_mu30e3.txt',
         '../pos_extremities/coupled_FoFo_375_CCEpara.txt',
         '../pos_extremities/bainite_FoFo_375.txt']

labels = [r"$\alpha'-\theta$" + u' ortho',
          # r'$\mu_C =$' + u' 16.0 kJ/mol',
          r'$\mu_C =$' + u' 20.0 kJ/mol',
          r'$\mu_C =$' + u' 23.2 kJ/mol',
          # r'$\mu_C =$' + u' 25.0 kJ/mol',
          r'$\mu_C =$' + u' 30.0 kJ/mol',
          r"$\alpha'-\theta$" + u' para',
          'no martensite']

for fname, label in zip(files, labels):
    phase_fractions = calculate_phase_fractions(fname)
    ax1.plot(phase_fractions.t, 100*phase_fractions.fer, label=label, lw=1)

ax1.legend(loc='lower right', fancybox=False, fontsize=10)
ax1.set_xlim(-1, 301)
ax1.set_ylim(-.5)
add_label(ax1, 'b)', px=.18, py=0)
ax1.set_xlabel('Isothermal holding (partitioning) time (s)')
ax1.set_ylabel(r'$f^{\alpha_b}$ (vol.%)')

##########################################################

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

y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
         Si=5.02504411E-2, Fe=9.4414085022e-1)
c0 = x2wp(3.34414e-02, y=y)
WBs = 1.5902

plot_cavg(files, labels=labels, ax=ax2, lw=1, y=y)
ax2.set_xlim(-1, 301)
ax2.set_ylim(-.01, 1.21)
ax2.text(298, c0, r'c$_0$', ha='right', va='bottom', size=10)
ax2.axhline(c0, ls=':', color='k', lw=1)
# ax2.axhline(WBs, ls=':', color='k', lw=1)
add_label(ax2, 'a)', px=.18, py=0)
ax2.legend(loc='lower right', fancybox=False, fontsize=10)

fig.tight_layout()
fout = 'cpartition_properties.svg'
fig.savefig(fout, bbox_inches='tight')
os.system('svg2pdf {}'.format(fout))

plt.show()
