# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def get_param_from_log(fname, param):
    idx, y = [], []
    with open(fname, 'r') as f:
        for line in f:
            if param in line:
                try:
                    line = line.strip().split(':')
                    idx.append(int(line[0]))
                except:
                    pass
                else:
                    args = line[1].strip().split(', ')
                    args = dict(arg.split('=') for arg in args)
                    s = args[param].strip('s')

                    y.append(float(s))

    return np.array(idx), np.array(y)


if __name__ == '__main__':
    import glob

    for fname in glob.glob('*.log'):
        idx, cavg = get_param_from_log(fname, 'cavg*')
        idx, t = get_param_from_log(fname, 't')

        plt.plot(t, cavg, label=fname)

    plt.legend()
    plt.show()
