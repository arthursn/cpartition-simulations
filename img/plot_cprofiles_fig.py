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

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax0, ax1, ax2, ax3 = axes.ravel()

    plt.subplots_adjust(wspace=.2, hspace=.35)

    ###########################################################

    # y, c0, and WBs for the high manganese, high carbon region (cell boundary)
    y = dict(Si=0.034356055114394574, Mn=0.003731497480722358,
             Cu=0.0027562013887263755)
    y['Fe'] = 1. - sum(y.values())
    c0 = x2wp(3.8553489675740495e-2, y=y)
    WBs = 1.5424

    ###########################################################

    labels = [('aus1', r'$\gamma$', 1),
              ('aus2', r'$\gamma$', 1),
              ('fer1', r'$\alpha_{b}$', 1),
              ('fer2', r'$\alpha_{b}$', 1),
              ('fer3', r'$\alpha_{b}$', 1)]

    p0 = CProfiles('bainite_FoFo_375', '../C_profiles')

    ax0.axhline(c0, ls=':', color='k', lw=1)
    ax0.axhline(WBs, ls=':', color='k', lw=1)
    ax0.text(-1.1, c0, r'c$_0$', ha='left', va='bottom', size=10)
    ax0.text(-1.15, WBs, 'WBs', ha='left', va='bottom', size=10)
    p0.plot_cprofiles(tlist=[1, 10, 100, 300], ax=ax0, func=lambda x: x2wp(x, y=y),
                      mirror=True, lw=1)

    ax0.set_xlim(-1.16, 1.16)
    ax0.set_ylim(-.02, 1.9)
    ax0.set_title('Bainite (cell boundary)', y=1.06, size=12)
    add_label(ax0, 'a)', py=0)
    ax0.set_xlabel(u'Position (μm)')
    ax0.set_ylabel('Carbon content (wt.%)')
    ax0.legend(loc='lower right', fancybox=False, fontsize=10)

    p0.label_phases(ax=ax0, t=300,
                    labels=labels,
                    mirror=True, size=12)

    ###########################################################

    # y, c0, and WBs for the average composition of the alloy
    y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
             Si=5.02504411E-2, Fe=9.4414085022e-1)
    c0 = x2wp(3.34414e-02, y=y)
    WBs = 1.5902

    ###########################################################

    labels = [('aus1', r'$\gamma$', 1),
              ('aus2', r'$\gamma$', 1),
              ('aust', r'$\gamma$', 1),
              ('fer1', r'$\alpha_{b}$', 1),
              ('fer2', r'$\alpha_{b}$', 1),
              ('mart', r"$\alpha' + \theta$", 1)]

    p1 = CProfiles('coupled_FoFo_375_CCEortho', '../C_profiles')

    ax1.axhline(c0, ls=':', color='k', lw=1)
    ax1.axhline(WBs, ls=':', color='k', lw=1)
    ax1.text(-1.1, c0, r'c$_0$', ha='left', va='bottom', size=10)
    ax1.text(-1.15, WBs, 'WBs', ha='left', va='bottom', size=10)
    p1.plot_cprofiles(tlist=[1, 10, 100], ax=ax1, func=lambda x: x2wp(x, y=y),
                      mirror=True, lw=1)

    ax1.set_xlim(-1.16, 1.16)
    ax1.set_ylim(-.02, 1.9)
    ax1.set_title(r"Martensite + bainite: $\alpha'-\theta$" + ' ortho', y=1.06, size=12)
    add_label(ax1, 'b)', py=0)
    ax1.set_xlabel(u'Position (μm)')
    ax1.set_ylabel('Carbon content (wt.%)')
    ax1.legend(loc='lower right', fancybox=False, fontsize=10)

    p1.label_phases(ax=ax1, t=100,
                    labels=labels,
                    mirror=True, size=12)

    ###########################################################

    p2 = CProfiles('coupled_FoFo_375_CCEpara', '../C_profiles')

    ax2.axhline(c0, ls=':', color='k', lw=1)
    ax2.axhline(WBs, ls=':', color='k', lw=1)
    ax2.text(-1.1, c0, r'c$_0$', ha='left', va='bottom', size=10)
    ax2.text(-1.15, WBs, 'WBs', ha='left', va='bottom', size=10)
    p2.plot_cprofiles(tlist=[1, 10, 100], ax=ax2, func=lambda x: x2wp(x, y=y),
                      mirror=True, lw=1)

    ax2.set_xlim(-1.16, 1.16)
    ax2.set_ylim(-.05, 4.3)
    ax2.set_title(r"Martensite + bainite: $\alpha'-\theta$" + ' para', y=1.06, size=12)
    add_label(ax2, 'c)', py=0)
    ax2.set_xlabel(u'Position (μm)')
    ax2.set_ylabel('Carbon content (wt.%)')
    ax2.legend(loc='lower right', fancybox=False, fontsize=10)

    p2.label_phases(ax=ax2, t=100,
                    labels=labels,
                    mirror=True, size=12)

    ###########################################################

    p3 = CProfiles('coupled_FoFo_375_mu23e3', '../C_profiles')

    ax3.axhline(c0, ls=':', color='k', lw=1)
    ax3.axhline(WBs, ls=':', color='k', lw=1)
    ax3.text(-1.1, c0, r'c$_0$', ha='left', va='bottom', size=10)
    ax3.text(-1.15, WBs, 'WBs', ha='left', va='bottom', size=10)
    p3.plot_cprofiles(tlist=[1, 10, 100], ax=ax3, func=lambda x: x2wp(x, y=y),
                      mirror=True, lw=1)

    ax3.set_xlim(-1.16, 1.16)
    ax3.set_ylim(-.02, 1.9)
    ax3.set_title(r"Martensite + bainite: $\mu_C=$" + '23.2 kJ/mol' +
                  r" ($c^{\gamma/\alpha'+\theta}$ = WBs)", y=1.06, size=12)
    add_label(ax3, 'd)', py=0)
    ax3.set_xlabel(u'Position (μm)')
    ax3.set_ylabel('Carbon content (wt.%)')
    ax3.legend(loc='lower right', fancybox=False, fontsize=10)

    p3.label_phases(ax=ax3, t=100,
                    labels=labels,
                    mirror=True, size=12)

    ###########################################################

    fout = 'carbon_profiles.svg'
    fig.savefig(fout, bbox_inches='tight')
    os.system('svg2pdf {}'.format(fout))

    plt.show()
    plt.close('all')
