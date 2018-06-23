import numpy as np
import matplotlib.pyplot as plt
from cpartition import FCC, BCC, Interface
from scipy.interpolate import UnivariateSpline

T_C = 375.

tdata_fcc = 'thermo/FoFo/TCFE8/375-fcc.txt'
tdata_bcc = 'thermo/FoFo/TCFE8/375-bcc.txt'

mart = BCC(T_C=T_C, tdata=tdata_bcc)
aust = FCC(T_C=T_C, tdata=tdata_fcc)
intf = Interface(domain1=mart, domain2=aust, type_int='fixed.fluxes')

# ph = mart
ph = aust
x, y = ph.chempot['X(C)'], ph.chempot['MU(C)']

_f = UnivariateSpline(np.log(x), y)
_f_prime = _f.derivative()


def f(x): return _f(np.log(x))


def f_prime(x): return _f_prime(np.log(x))/x


def M(x):
    return ph.D(x)/(x*f_prime(x))


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(x, ph.D(x), 'r-')
ax2.plot(x, M(x), 'b-')

ax1.set_xlabel('Carbon content (at.%)')
ax1.set_ylabel(r'Carbon diffusivity ($\mu m^2/s$)')
ax2.set_ylabel(r'Carbon mobility (??)')

plt.show()
