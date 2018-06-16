#!/usr/bin/env python

import numpy as np
import time

import os
import sys
from cpartition import *
from scipy.interpolate import interp1d

basename = os.path.basename(__file__).replace('.py', '')

c0 = 3.34414e-02
T_C = 375.

dt = 5e-4
total_time = 1
n_time = int(total_time/dt)

tdata_fcc = 'thermo/FoFo/TCFE8/375-fcc.txt'
tdata_bcc = 'thermo/FoFo/TCFE8/375-bcc.txt'

mart = BCC(T_C=T_C, dt=dt, z=np.linspace(-1.16, -.66, 50), c0=c0,
           n_time=n_time, tdata=tdata_bcc)
aust = FCC(T_C=T_C, dt=dt, z=np.linspace(-.66, 0, 400), c0=c0,
           n_time=n_time, tdata=tdata_fcc)

intf = Interface(domain1=mart, domain2=aust, type_int='fixed.flux')

log = SimulationLog(basename)
log.set_domains([('mart', mart), ('aust', aust)])
log.set_interfaces([('intf', intf)])
log.set_conditions(c0, T_C, total_time, n_time)
log.initialize(1, False)

for i in range(n_time):
    intf.comp(poly_deg=2)
    mart.FDM_implicit(bcn=(1., 0, 0, intf.ci_bcc))
    aust.FDM_implicit(bc0=(1., 0, 0, intf.ci_fcc))

    mart.update_grid()
    aust.update_grid()

    log.print(i)

log.close()

log.save_cprofiles()
log.save_properties('cavg')
log.save_properties('ci*')
log.save_properties('s*')

##############

dt = 5e-2
total_time = 1000
n_time = int(total_time/dt)

z = np.linspace(mart.z[0], mart.z[-1], 10)
c = interp1d(mart.z, mart.c)(z)
mart = BCC(T_C=T_C, dt=dt, z=z, c=c,
           n_time=n_time, tdata=tdata_bcc)

z = np.linspace(aust.z[0], aust.z[-1], 50)
c = interp1d(aust.z, aust.c)(z)
aust = FCC(T_C=T_C, dt=dt, z=z, c=c,
           n_time=n_time, tdata=tdata_fcc)


intf = Interface(domain1=mart, domain2=aust, type_int='fixed.flux')

log = SimulationLog(basename + '_resume')
log.set_domains([('mart', mart), ('aust', aust)])
log.set_interfaces([('intf', intf)])
log.set_conditions(c0, T_C, total_time, n_time)
log.initialize(100, False)

for i in range(n_time):
    intf.comp(poly_deg=2)
    mart.FDM_implicit(bcn=(1., 0, 0, intf.ci_bcc))
    aust.FDM_implicit(bc0=(1., 0, 0, intf.ci_fcc))

    mart.update_grid()
    aust.update_grid()

    log.print(i)


log.close()

log.save_cprofiles()
log.save_properties('cavg')
log.save_properties('ci*')
log.save_properties('s*')
