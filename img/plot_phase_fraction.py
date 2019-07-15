import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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


fig, ax = plt.subplots(figsize=(5, 4))


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
          'austempering']

for fname, label in zip(files, labels):
    phase_fractions = calculate_phase_fractions(fname)
    ax.plot(phase_fractions.t, phase_fractions.fer, label=label)

ax.legend(fancybox=False)
ax.set_xlim(-.5, 200.5)

plt.show()


# coupled = calculate_phase_fractions('../pos_extremities/coupled_FoFo_375_mu30e3.txt')
# bainite = calculate_phase_fractions('../pos_extremities/bainite_FoFo_375.txt')

# fig1, ax1 = plt.subplots(figsize=(5, 4))
# ax1.plot(coupled.t, coupled.fer, label='coupled')
# ax1.plot(bainite.t, bainite.fer, label='bainite')
# ax1.legend(fancybox=False)

# ax1.set_xlim(0, 200)
# ax1.set_ylim(0)
# ax1.set_xlabel('Time (s)')
# ax1.set_ylabel(r'$f^{\alpha_b}$')
# fig1.tight_layout()

# fig2, ax2 = plt.subplots(figsize=(5, 4))
# ax2.plot(coupled.t, coupled.fer/coupled.aus[0], label='coupled')
# ax2.plot(bainite.t, bainite.fer/bainite.aus[0], label='bainite')
# ax2.legend(fancybox=False)

# ax2.set_xlim(0, 200)
# ax2.set_ylim(0, 1)
# ax2.set_xlabel('Time (s)')
# ax2.set_ylabel(r'$f^{\alpha_b}/f^\gamma_0$')

# ax21 = ax2.twinx()
# ax21.set_ylim(1, 0)
# ax21.set_ylabel(r'$f^{\gamma}/f^\gamma_0$')

# fig2.tight_layout()

# plt.show()
