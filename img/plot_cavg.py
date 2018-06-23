# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle
from cpartition import x2wp


def x2vcem(x):
    vmbcc = 6.75867413e-6
    vmcem = 5.7851763e-6
    return x*vmcem/(x*vmcem + (1 - x)*vmbcc)


def plot_cavg(files, labels=None, styles=None, ax=None,
              plotcarbides=False, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    if labels:
        cylbls = cycle(labels)

    if styles:
        cysty = cycle(styles)

    x = kwargs.pop('x', None)
    w = kwargs.pop('w', None)
    y = kwargs.pop('y', None)

    for fname in files:
        df = pd.read_table(fname, sep=' ')

        args = []
        if labels:
            kwargs.update(dict(label=next(cylbls)))
        if styles:
            args.append(next(cysty))

        ax.plot(df.t, x2wp(df['mart.cavg'], x=x, w=w, y=y), *args, **kwargs)

    ax.set_xlim(-.5, 100.5)

    yrng = np.array([-1e-4, 4.5e-2])
    ax.set_ylim(x2wp(yrng, x=x, w=w, y=y))

    ax.set_xlabel('Partitioning time (s)')
    ax.set_ylabel(r"Carbon in $\alpha' + \theta$ (wt.%)")

    if plotcarbides:
        ax2 = ax.twinx()
        ax2.set_ylim(100*x2vcem(yrng/.25))
        ax2.set_ylabel("Carbides fraction (vol.%)")
        ax.legend(loc='lower right', fancybox=False)

        return ax, ax2
    else:
        return ax


if __name__ == '__main__':
    import sys
    from matplotlib import rcParams

    rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': 'Arial',
                     'font.size': 13})

    args = sys.argv[1:]
    if len(args) > 0:
        plot_cavg(files=args, labels=args)
        plt.show()
    else:
        print('Nothing to plot')
