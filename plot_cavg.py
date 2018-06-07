# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

from matplotlib import rcParams
import matplotlib.pyplot as plt

from itertools import cycle
from periodictable import elements


def x2w(xC, y=dict(Cu=3.55354266E-3, Mn=2.05516602E-3, Si=5.02504411E-2)):
    x = {k: v*(1-xC) for k, v in y.items()}
    x['C'] = xC
    if 'Fe' not in x.keys():
        x['Fe'] = 1 - sum(x.values())

    w = {k: v*elements.symbol(k).mass for k, v in x.items()}

    return w['C']/sum(w.values())


def x2vcem(x):
    vmbcc = 6.75867413e-6
    vmcem = 5.7851763e-6
    return x*vmcem/(x*vmcem + (1 - x)*vmbcc)


def plot_cavg(files, labels=None, styles=None, ax=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    if labels:
        cylbls = cycle(labels)

    if styles:
        cysty = cycle(styles)

    for fname in files:
        df = pd.read_table(fname, sep=' ')

        args = []
        if labels:
            kwargs.update(dict(label=next(cylbls)))
        if styles:
            args.append(next(cysty))

        ax.plot(df.t, 100*x2w(df.mart), *args, **kwargs)

    ax2 = ax.twinx()
    ax.set_xlim(-.5, 100.5)

    yrng = np.array([-1e-4, 4.5e-2])
    ax.set_ylim(100*x2w(yrng))
    ax2.set_ylim(100*x2vcem(yrng/.25))

    ax.set_xlabel('Partitioning time (s)')
    ax.set_ylabel(r"Carbon in $\alpha' + \theta$ (wt.%)")
    ax2.set_ylabel("Carbides fraction (vol.%)")
    ax.legend(loc='lower right', fancybox=False)

    return ax


if __name__ == '__main__':
    rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': 'Arial',
                     'font.size': 13,
                     'mathtext.fontset': 'stix'})

    import sys

    if len(sys.argv) > 1:
        files = sys.argv[1:]

        plot_cavg(files=files, labels=files)
        plt.show()
    else:
        print('Nothing to plot')
