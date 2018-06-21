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

control_itsteps = ControlIterationSteps([5e-5, 5e-4, 5e-3, 5e-2], [0, .2, 2, 20, 1000])
total_time = control_itsteps.total_time
n_time = control_itsteps.ntime
dt = control_itsteps.dt
each = 200
control_itsteps.print_summary()

tdata_fcc = 'thermo/FoFo/TCFE8/375-fcc.txt'
tdata_bcc = 'thermo/FoFo/TCFE8/375-bcc.txt'

mart = BCC(T_C=T_C, dt=dt, z=np.linspace(-1.16, -.66, 20), c0=c0,
           n_time=n_time, tdata=tdata_bcc,
           type_D='carbides', cmax_bcc=5.4e-4, c_carbide=.25)
aust = FCC(T_C=T_C, dt=dt, z=np.linspace(-.66, 0, 100), c0=c0,
           n_time=n_time, tdata=tdata_fcc)

intf = Interface(domain1=mart, domain2=aust, type_int='fixed.flux')

# fixed composition set by CCEtheta in austenite at the interface
muC = 35511.1
cCCEtheta = aust.mu2x['C'](muC)
intf.ci_fcc = cCCEtheta
print('muC={:}, ci_fcc={:}\n'.format(muC, intf.ci_fcc))

log = SimulationLog(basename)
log.set_domains([('mart', mart), ('aust', aust)])
log.set_interfaces([('intf', intf)])
log.set_conditions(c0, T_C, total_time_1, n_time)
log.initialize(False)

for i in control_itsteps.itlist:
    if i in control_itsteps.itstepi and i > 0:
        control_itsteps.next_itstep()
        dt = control_itsteps.dt

        mart.dt = dt
        aust.dt = dt

    J, = intf.flux('fcc')

    grad = -J*mart.dz/mart.D(mart.c[-1])
    if mart.FDM_implicit(bcn=(1.5, -2., .5, grad), lowerbound=5.4e-4):
        intf.comp(poly_deg=3)
        mart.FDM_implicit(bcn=(1., 0, 0, intf.ci_bcc))
        aust.FDM_implicit(bc0=(1., 0, 0, intf.ci_fcc))
    else:
        aust.FDM_implicit(bc0=(1., 0, 0, cCCEtheta))

    mart.update_grid(i)
    aust.update_grid(i)

    log.printit(i, each)

# log.close()

# dt = 5e-3
# total_time_2 = 1000
# n_time = int((total_time_2 - total_time_1)/dt)
# each = 200

# z = np.linspace(mart.z[0], mart.z[-1], 10)
# c = interp1d(mart.z, mart.c)(z)
# mart = BCC(T_C=T_C, dt=dt, z=z, c=c,
#            t0=total_time_1, tdata=tdata_bcc)

# z = np.linspace(aust.z[0], aust.z[-1], 100)
# c = interp1d(aust.z, aust.c)(z)
# aust = FCC(T_C=T_C, dt=dt, z=z, c=c,
#            t0=total_time_1, tdata=tdata_fcc)


# intf = Interface(domain1=mart, domain2=aust, type_int='fixed.flux')

# log.set_domains([('mart', mart), ('aust', aust)])
# log.set_interfaces([('intf', intf)])
# log.set_conditions(c0, T_C, total_time_1, n_time, reset=False)
# log.initialize(False, mode='a')

# def new_crit(it, each):
#     return (it - it1 + int(total_time_1/dt)) % each == 0

# for it2 in range(it1 + 1, it1 + n_time + 1):
#     intf.comp(poly_deg=2)
#     mart.FDM_implicit(bcn=(1., 0, 0, intf.ci_bcc))
#     aust.FDM_implicit(bc0=(1., 0, 0, intf.ci_fcc))

#     mart.update_grid(it2)
#     aust.update_grid(it2)

#     log.printit(it2, each, criteria=new_crit)


log.close()

log.save_cprofiles()
log.save_properties('cavg')
log.save_properties('ci*')
log.save_properties('s*')
