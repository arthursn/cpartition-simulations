# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt

import os
from cpartition import x2wp, label, CProfiles

from plot_cavg import plot_cavg


class PlotParams(object):
    def __init__(self, basename):
        self.basename = basename
        self.cprofiles = None

    @property
    def outappend(self):
        outappend = ''
        if self.log is True:
            outappend += '_log'
        if self.tracking is True:
            outappend += '_tracking'
        return outappend

    def load_file(self):
        try:
            fname = self.basename + '_profiles.txt'
            self.cprofiles = CProfiles(fname)
            self.cprofiles.load_file()
        except:
            print('File "{}" does not exist'.format(fname))
            raise
        else:
            print('File "{}" successfully loaded'.format(fname))

    def plot_profiles(self, t_set, ax=None, func=lambda x: x, **kwargs):
        if not ax:
            fig, ax = plt.subplots(figsize=self.figsize)
        else:
            fig = ax.get_figure()

        if not self.cprofiles:
            self.load_file()

        if self.cprofiles:
            mirror = kwargs.pop('mirror', False)
            xlim = kwargs.pop('xlim', None)
            ylim = kwargs.pop('ylim', None)
            log = kwargs.pop('log', False)
            loc = kwargs.pop('loc', 'upper left')

            for t in t_set:
                idx = int(t/self.cprofiles.tstep) - 1

                if idx > 0:
                    z, c = self.cprofiles.zz[idx], self.cprofiles.cc[idx]

                    if mirror:
                        z = np.hstack([z, 2*z[-1] - z[::-1]])
                        c = np.hstack([c, c[::-1]])
                    ax.plot(z, func(c), label=label(t), **kwargs)

                    try:
                        if log:
                            ax.set_yscale('log')
                    except:
                        pass

                    ax.set_xlim(xlim)
                    ax.set_ylim(ylim)
                    ax.set_xlabel(u'Position (µm)')
                    ax.set_ylabel('Carbon content (wt.%)')
                    ax.legend(loc=loc, fancybox=False)
                else:
                    print('t = {} < tmin'.format(t))


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
    p1 = PlotParams('coupled_FoFo_375_CCEortho')
    ax1.axhline(c0, ls=':', color='k', lw=1)
    ax1.axhline(WBs, ls=':', color='k', lw=1)
    p1.plot_profiles(t_set=[1, 10, 50], ax=ax1, func=lambda x: x2wp(x, y=y),
                     mirror=True, xlim=(-1.16, 1.16), ylim=(-.02, 1.9), lw=1)
    add_label(ax1, 'a)', py=0)
    ax1.text(.98, .98, r"$\alpha'-\theta$" + ' ortho',
             ha='right', va='top', transform=ax1.transAxes)

    ############
    p2 = PlotParams('coupled_FoFo_375_CCEpara')
    ax2.axhline(c0, ls=':', color='k', lw=1)
    ax2.axhline(WBs, ls=':', color='k', lw=1)
    p2.plot_profiles(t_set=[1, 10, 100], ax=ax2, func=lambda x: x2wp(x, y=y),
                     mirror=True, xlim=(-1.16, 1.16), ylim=(-.05, 4.1), lw=1)
    add_label(ax2, 'b)', py=0)
    ax2.text(.98, .98, r"$\alpha'-\theta$" + ' para',
             ha='right', va='top', transform=ax2.transAxes)

    ############
    p3 = PlotParams('coupled_FoFo_375_mu23e3')
    ax3.axhline(c0, ls=':', color='k', lw=1)
    ax3.axhline(WBs, ls=':', color='k', lw=1)
    p3.plot_profiles(t_set=[1, 10, 100], ax=ax3, func=lambda x: x2wp(x, y=y),
                     mirror=True, xlim=(-1.16, 1.16), ylim=(-.02, 1.9), lw=1)
    add_label(ax3, 'c)', py=0)
    ax3.text(.98, .98, r"$\mu_C=$" + '23.2 kJ/mol' + r' ($c^\gamma_{int}$ = WBs)',
             ha='right', va='top', transform=ax3.transAxes)

    ############
    files = ['C_avg/coupled_FoFo_375_CCEortho.txt',
             'C_avg/coupled_FoFo_375_mu20e3.txt',
             'C_avg/coupled_FoFo_375_mu23e3.txt',
             'C_avg/coupled_FoFo_375_mu30e3.txt',
             'C_avg/coupled_FoFo_375_CCEpara.txt']
    labels = [r"$\alpha'-\theta$" + u' ortho',
              r'$\mu_C =$' + u' 20.0 kJ/mol',
              r'$\mu_C =$' + u' 23.2 kJ/mol',
              r'$\mu_C =$' + u' 30.0 kJ/mol',
              r"$\alpha'-\theta$" + u' para']
    # styles = ['k--', 'k-.', 'k-']

    plot_cavg(files, labels, ax=ax4, lw=1)
    ax4.axhline(c0, ls=':', color='k', lw=1)
    ax4.axhline(WBs, ls=':', color='k', lw=1)
    add_label(ax4, 'd)', py=0)

    fout = 'img/cpartition.pdf'
    fig.savefig(fout, bbox_inches='tight')

    plt.show()
    plt.close('all')
