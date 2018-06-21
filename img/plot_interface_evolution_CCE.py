# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import pandas as pd
from cpartition import FCC, BCC, Interface, WBs, x2wp

import matplotlib
matplotlib.rc('font', **{'family':'sans-serif', 'sans-serif':['Arial'], 'size': 13})

# Site fraction of substitutional elements in cast iron's austenite
T_C = 375
y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
         Si=5.02504411E-2, Fe=9.4414085022e-1)

cint = pd.read_table('../C_extremities/coupled_FoFo_375_CCE.txt', sep=' ')
pos = pd.read_table('../pos_extremities/coupled_FoFo_375_CCE.txt', sep=' ')

# For calculating WBs
aust = FCC(T_C=T_C, tdata='../thermo/FoFo/TCFE8/375-fcc.txt')
ferr = BCC(T_C=T_C, tdata='../thermo/FoFo/TCFE8/375-bcc.txt', E=WBs(T_C))
intf = Interface(domain1=aust, domain2=ferr, type_int='mobile.eq')
# WBs composition
_, cwbs = intf.comp()
cwbs = x2wp(cwbs, y=y)

# Driving force
DF = intf.chem_driving_force(ci_bcc=cint['fer1.ci0'], ci_fcc=cint['aus1.cin'])

pos0 = pos['aus1.sn'].values[0]
tf = pos['t'].values[-1]

# Ploting interface position and composition

# Setting pre-plot parameters
fig1, ax1 = plt.subplots(figsize=(4, 4))

# composition
ax1.plot(cint['t'], x2wp(cint['aus1.cin'], y=y),
         'k-', label=r'$c^\gamma_{int_1}$')
ax1.axhline(cwbs, c='k', ls=':')

ax1.text(x=990, y=cwbs - .01, s=u'WBs ({:.2f} %)'.format(cwbs),
         fontsize=16, ha='right', va='top')
ax1.set_ylim(cwbs - .17, cwbs + .17)
ax1.set_xlabel(u'Tempo (s)')
ax1.set_ylabel(r'$c^\gamma_{int_1}$ (% massa)')

# position
ax2 = ax1.twinx()
ax2.plot(pos['t'], pos['aus1.sn'], 'k--', label=u'Posição da interface')

ax2.set_xlim(1e0, tf)
ax2.set_xscale('log')
ax2.set_ylim(pos0 - .1, pos0 + .25)
ax2.set_ylabel(u'Posição da interface (μm)')
ax2.invert_yaxis()

# legend
ax1.legend(loc='lower left', fancybox=False)
ax2.legend(loc='upper right', fancybox=False)

# Ploting driving force

fig3, ax3 = plt.subplots(figsize=(4, 4))

ax3.plot(cint['t'], DF, 'k-', label="Força motriz")
ax3.set_xlim(1e0, tf)
ax3.set_xscale('log')
ax3.set_ylim(-20, 20)
ax3.set_xlabel('Tempo (s)')
ax3.set_ylabel(u'Força motriz (J/mol)')
ax3.legend(loc='upper right', fancybox=False)
ax3.axhline(0, c='k', ls=':')

f1 = '/home/arthur/tese/img/cpartition/aus1fer1_interface_comp.svg'
f2 = '/home/arthur/tese/img/cpartition/aus1fer1_interface_DF.svg'

fig1.savefig(f1, bbox_inches='tight')
fig3.savefig(f2, bbox_inches='tight')

import os
os.system('svg2pdf ' + f1)
os.system('svg2pdf ' + f2)

plt.show()
plt.close()
