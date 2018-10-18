# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt

import os
from cpartition import x2wp, label, CProfiles

# script in the local folder
from plot_cavg import plot_cavg


def add_label(ax, label, px=.15, py=.1, size=20, **kwargs):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x = float(min(xlim) - np.diff(xlim)*px)
    y = float(max(ylim) + np.diff(ylim)*py)
    ax.text(x=x, y=y, s=label, size=size, **kwargs)


if __name__ == '__main__':
    rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': 'Arial',
                     'font.size': 13,
                     'mathtext.fontset': 'stix'})

    y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
             Si=5.02504411E-2, Fe=9.4414085022e-1)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax1, ax2, ax3, ax4 = axes.ravel()

    plt.subplots_adjust(wspace=.3, hspace=.3)

    c0 = x2wp(3.34414e-02, y=y)
    WBs = 1.5902

    ############
    labels = [('aus1', r'$\gamma$', 1),
              ('aus2', r'$\gamma$', 1),
              ('aust', r'$\gamma$', 1),
              ('fer1', r'$\alpha_{b}$', 1),
              ('fer2', r'$\alpha_{b}$', 1),
              ('mart', r"$\alpha' + \theta$", 1)]

    p1 = CProfiles('coupled_FoFo_375_CCEortho', '../C_profiles')

    ax1.axhline(c0, ls=':', color='k', lw=1)
    ax1.axhline(WBs, ls=':', color='k', lw=1)
    p1.plot_cprofiles(tlist=[1, 10, 100], ax=ax1, func=lambda x: x2wp(x, y=y),
                      mirror=True, lw=1)

    ax1.set_xlim(-1.16, 1.16)
    ax1.set_ylim(-.02, 1.9)
    ax1.text(.95, .95, r"$\alpha'-\theta$" + ' ortho',
             ha='right', va='top', size=10,
             transform=ax1.transAxes,
             backgroundcolor='white')
    add_label(ax1, 'a)', py=0)
    ax1.set_xlabel(u'Position (μm)')
    ax1.set_ylabel('Carbon content (wt.%)')
    ax1.legend(loc='lower right', fancybox=False, fontsize=10)

    p1.label_phases(ax=ax1, t=1,
                    labels=labels,
                    mirror=True, size=12)

    ############
    p2 = CProfiles('coupled_FoFo_375_CCEpara', '../C_profiles')

    ax2.axhline(c0, ls=':', color='k', lw=1)
    ax2.axhline(WBs, ls=':', color='k', lw=1)
    p2.plot_cprofiles(tlist=[1, 10, 100], ax=ax2, func=lambda x: x2wp(x, y=y),
                      mirror=True, lw=1)

    ax2.set_xlim(-1.16, 1.16)
    ax2.set_ylim(-.05, 4.3)
    ax2.text(.95, .95, r"$\alpha'-\theta$" + ' para',
             ha='right', va='top', size=10,
             transform=ax2.transAxes,
             backgroundcolor='white')
    add_label(ax2, 'b)', py=0)
    ax2.set_xlabel(u'Position (μm)')
    ax2.set_ylabel('Carbon content (wt.%)')
    ax2.legend(loc='lower right', fancybox=False, fontsize=10)

    p2.label_phases(ax=ax2, t=1,
                    labels=labels,
                    mirror=True, size=12)

    ############
    p3 = CProfiles('coupled_FoFo_375_mu23e3', '../C_profiles')

    ax3.axhline(c0, ls=':', color='k', lw=1)
    ax3.axhline(WBs, ls=':', color='k', lw=1)
    p3.plot_cprofiles(tlist=[1, 10, 100], ax=ax3, func=lambda x: x2wp(x, y=y),
                      mirror=True, lw=1)

    ax3.set_xlim(-1.16, 1.16)
    ax3.set_ylim(-.02, 1.9)
    ax3.text(.95, .95, r"$\mu_C=$" + '23.2 kJ/mol' + r' ($c^\gamma_{int}$ = WBs)',
             ha='right', va='top', size=10,
             transform=ax3.transAxes,
             backgroundcolor='white')
    add_label(ax3, 'c)', py=0)
    ax3.set_xlabel(u'Position (μm)')
    ax3.set_ylabel('Carbon content (wt.%)')
    ax3.legend(loc='lower right', fancybox=False, fontsize=10)

    p3.label_phases(ax=ax3, t=1,
                    labels=labels,
                    mirror=True, size=12)

    # Plot average C composition

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

    plot_cavg(files, labels, ax=ax4, lw=1, y=y)
    ax4.axhline(c0, ls=':', color='k', lw=1)
    ax4.axhline(WBs, ls=':', color='k', lw=1)
    add_label(ax4, 'd)', py=0)
    ax4.legend(loc='lower right', fancybox=False, fontsize=10)

    fout = 'cpartition.svg'
    fig.savefig(fout, bbox_inches='tight')
    os.system('svg2pdf {}'.format(fout))

    plt.show()
    plt.close('all')
