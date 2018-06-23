#!/usr/bin/env python

import numpy as np
import time

import os
import sys
from cpartition import *

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

mart = BCC(T_C=T_C, dt=dt, z=np.linspace(-1.16, -.66, 50), c0=c0,
           n_time=n_time, tdata=tdata_bcc,
           type_D='carbides', cmax_bcc=5.4e-4, c_carbide=.3)
aust = FCC(T_C=T_C, dt=dt, z=np.linspace(-.66, 0, 200), c0=c0,
           n_time=n_time, tdata=tdata_fcc)

intf = Interface(domain1=mart, domain2=aust, type_int='fixed.fluxes')

# fixed composition set by CCEtheta in austenite at the interface
muC = 30e3
cCCEtheta = aust.mu2x['C'](muC)
intf.ci_fcc = cCCEtheta
print('muC={:}, ci_fcc={:}\n'.format(muC, intf.ci_fcc))

log = SimulationLog(basename)
log.set_domains([('mart', mart), ('aust', aust)])
log.set_interfaces([('intf', intf)])
log.set_conditions(c0, T_C, total_time, n_time)
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

log.close()

log.save_cprofiles()
log.save_properties('cavg')
log.save_properties('ci*')
log.save_properties('s*')
