#!/usr/bin/env python

import numpy as np
import time

import os
import sys
newdir = '/home/arthur/Dropbox/python'
if newdir not in sys.path:
    sys.path.insert(1, newdir)
from cpartition import *
import matplotlib.pyplot as plt

"""Coupled model for Fe-0.76C alloy"""

c0 = 3.34414e-02
T_C = 375.

n_time = 2000
total_time = 1.
dt = total_time/n_time
t = (np.arange(n_time) + 1)*dt

tdata_fcc = 'thermo/FoFo/TCFE8/375-fcc.txt'
tdata_bcc = 'thermo/FoFo/TCFE8/375-bcc.txt'

mart = BCC(T_C=T_C, tdata=tdata_bcc)
aust = FCC(T_C=T_C, tdata=tdata_fcc)

intf = Interface(domain1=mart, domain2=aust, type_int='fixed.fluxes')

ph = aust
# ph = mart
x, y = ph.chempot['X(C)'], ph.chempot['MU(C)']
_f = UnivariateSpline(np.log(x), y)
_f_prime = _f.derivative()
f = lambda x: _f(np.log(x))
f_prime = lambda x: _f_prime(np.log(x))/x

def M(x):
    return ph.D(x)/(x*f_prime(x))

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(x, ph.D(x), 'r-')
ax2.plot(x, M(x), 'b-')
plt.show()
