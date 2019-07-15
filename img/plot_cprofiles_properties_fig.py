import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt

import os
from cpartition import x2wp, label, CProfiles

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


def calculate_phase_fractions(fname_position):
    df = pd.read_csv(fname_position, delim_whitespace=True)

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


def calculate_average_carbon(fname_position, fname_cavg):
    df = pd.read_csv(fname_position, delim_whitespace=True)

    domain_phase_map = dict()

    for colname in df.columns[1:]:
        domain = colname.split('.')[0]
        phase = ''.join([i for i in domain if not i.isdigit()])

        if phase not in domain_phase_map.keys():
            domain_phase_map[phase] = set([domain])
        else:
            domain_phase_map[phase].add(domain)

    df_cavg = pd.read_csv(fname_cavg, delim_whitespace=True)

    average_carbon = pd.DataFrame(0, index=df.index, columns=domain_phase_map.keys())

    for phase, domains in domain_phase_map.items():
        length_phase = 0.
        for domain in domains:
            length_domain = (df['{}.sn'.format(domain)] - df['{}.s0'.format(domain)]).abs()
            length_phase += length_domain
            average_carbon[phase] += df_cavg['{}.cavg'.format(domain)]*length_domain

        average_carbon[phase] = average_carbon[phase]/length_phase

    average_carbon['t'] = df['t']
    cols = average_carbon.columns.tolist()
    average_carbon = average_carbon[cols[-1:] + cols[:-1]]

    return average_carbon


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(5, 9))
# plt.subplots_adjust(hspace=.1)

##########################################################

files_cavg = ['../C_avg/coupled_FoFo_375_CCEortho.txt',
              '../C_avg/coupled_FoFo_375_mu20e3.txt',
              '../C_avg/coupled_FoFo_375_mu23e3.txt',
              '../C_avg/coupled_FoFo_375_mu30e3.txt',
              '../C_avg/coupled_FoFo_375_CCEpara.txt',
              '../C_avg/bainite_FoFo_375.txt']

files_position = ['../pos_extremities/coupled_FoFo_375_CCEortho.txt',
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

##########################################################

# c0 = x2wp(3.34414e-02, y=y)
# WBs = 1.5902

for fname_position, fname_cavg, label in zip(files_position, files_cavg, labels):
    average_carbon = calculate_average_carbon(fname_position, fname_cavg)

    if 'bainite' not in fname_position:
        y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
                 Si=5.02504411E-2, Fe=9.4414085022e-1)
        ax1.plot(average_carbon.t, x2wp(average_carbon.mart, y=y), label=label, lw=1)
    else:
        y = dict(Si=0.034356055114394574, Mn=0.003731497480722358,
                 Cu=0.0027562013887263755)
        y['Fe'] = 1. - sum(y.values())

    ax2.plot(average_carbon.t, x2wp(average_carbon.aus, y=y), label=label, lw=1)


ax1.set_xlabel('Isothermal holding (partitioning) time (s)')
ax1.set_ylabel(r"$\overline{c^{\alpha' + \theta}}$ (wt.%)")
ax1.set_xlim(-2, 602)
ax1.set_ylim(-.01, 1.31)
# ax1.text(298, c0, r'c$_0$', ha='right', va='bottom', size=10)
# ax1.axhline(c0, ls=':', color='k', lw=1)
add_label(ax1, 'a)', px=.2, py=0)

ax2.set_xlabel('Isothermal holding (partitioning) time (s)')
ax2.set_ylabel(r"$\overline{c^{\gamma}}$ (wt.%)")
ax2.set_xlim(-2, 602)
ax2.set_ylim(.7, 2.)
# ax2.text(298, c0, r'c$_0$', ha='right', va='bottom', size=10)
# ax2.axhline(c0, ls=':', color='k', lw=1)
add_label(ax2, 'b)', px=.2, py=0)

##########################################################

lines = []

for fname, label in zip(files_position, labels):
    phase_fractions = calculate_phase_fractions(fname)
    lines += ax3.plot(phase_fractions.t, 100*phase_fractions.fer, label=label, lw=1)

ax3.set_xlim(-2, 602)
ax3.set_ylim(-.5)
add_label(ax3, 'c)', px=.2, py=0)
ax3.set_xlabel('Isothermal holding (partitioning) time (s)')
ax3.set_ylabel(r'$f^{\alpha_b}$ (vol.%)')

##########################################################

ax1.legend(lines, labels, loc='lower right', fancybox=False, fontsize=10)

fig.tight_layout()
fout = 'cpartition_properties.svg'
fig.savefig(fout, bbox_inches='tight')
os.system('svg2pdf {}'.format(fout))

plt.show()
