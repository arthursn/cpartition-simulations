#!/usr/bin/env python

import numpy as np
import time

import sys
import os
sys.path.insert(1, '/home/arthur/Dropbox/python')
from cpartition import *

from scipy.interpolate import interp1d


T_C = 375.

tdata_fcc = 'thermo/FoFo/TCFE8/375-fcc.txt'
tdata_bcc = 'thermo/FoFo/TCFE8/375-bcc.txt'
# tdata_fcc = 'thermo/FoFo/TCFE0/375-FCC.TXT'
# tdata_bcc = 'thermo/FoFo/TCFE0/375-BCC.TXT'
aust = FCC(T_C=T_C, tdata=tdata_fcc)
ferr = BCC(T_C=T_C, tdata=tdata_bcc, E=WBs(T_C))

def g(x): return aust.f(x) - ferr.f(x)
lo, hi = max(min(aust.muC), min(ferr.muC)), min(max(aust.muC), max(ferr.muC))

muC = bisect(g, lo, hi, xtol=1e-3)
cferr = ferr.mu2x['C'](muC)
caust = aust.mu2x['C'](muC)

print('muC={:}, caust={:}\n'.format(muC, caust))

# import matplotlib.pyplot as plt
# plt.plot(aust.muZ, aust.muC)
# plt.plot(ferr.muZ, ferr.muC)
# plt.show()

