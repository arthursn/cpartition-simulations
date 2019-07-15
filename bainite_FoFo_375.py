import numpy as np
import time

import os
import sys
from cpartition import *

basename = os.path.basename(__file__).replace('.py', '')

# wc0 = 0.76e-2
# c0 = 3.34414e-02
# Carbon in the cell boundary
c0 = 3.8553489675740495e-2
T_C = 375.

control_itsteps = ControlIterationSteps([5e-5, 5e-4, 5e-3, 5e-2], [0, .2, 2, 20, 1000])
total_time = control_itsteps.total_time
n_time = control_itsteps.ntime
dt = control_itsteps.dt
each = 20
control_itsteps.print_summary()

tdata_fcc = 'thermo/FoFo/TCFE9/fofo_cell_boundary/375-FCC.TXT'
tdata_bcc = 'thermo/FoFo/TCFE9/fofo_cell_boundary/375-BCC.TXT'

fer2_pos = -1.16/3.
fer1 = BCC(T_C=T_C, dt=dt, z=np.linspace(-1.16, -1.16, 10), c0=0.,
           tdata=tdata_bcc, E=WBs(T_C))
aus1 = FCC(T_C=T_C, dt=dt, z=np.linspace(-1.16, fer2_pos, 100), c0=c0,
           tdata=tdata_fcc)
fer2 = BCC(T_C=T_C, dt=dt, z=np.linspace(fer2_pos, fer2_pos, 10), c0=0.,
           tdata=tdata_bcc, E=WBs(T_C))
aus2 = FCC(T_C=T_C, dt=dt, z=np.linspace(fer2_pos, 0, 100), c0=c0,
           tdata=tdata_fcc)

int1 = Interface(domain1=fer1, domain2=aus1, type_int='mobile.mmode')
int2 = Interface(domain1=aus1, domain2=fer2, type_int='mobile.mmode')
int3 = Interface(domain1=fer2, domain2=aus2, type_int='mobile.mmode')

fer1.c[:] = int1.CCE(c0)
fer2.c[:] = int2.CCE(c0)

log = SimulationLog(basename)
log.set_domains([('fer1', fer1), ('aus1', aus1),
                 ('fer2', fer2), ('aus2', aus2)])
log.set_interfaces([('int1', int1), ('int2', int2),
                    ('int3', int3)])
log.set_conditions(c0, T_C, total_time, n_time)
log.initialize(False)

upgrade = False

for i in control_itsteps.itlist:
    if i in control_itsteps.itstepi and i > 0:
        control_itsteps.next_itstep()
        dt = control_itsteps.dt

        fer1.dt = dt
        aus1.dt = dt
        fer2.dt = dt
        aus2.dt = dt

    # interface velocities at the mobile interfaces
    int1.v = 1e6*int1.chem_driving_force()*int1.M()/fer1.Vm
    int1.comp(poly_deg=3)
    int2.v = 1e6*int2.chem_driving_force()*int2.M()/fer2.Vm
    int2.comp(poly_deg=3)
    int3.v = 1e6*int3.chem_driving_force()*int3.M()/fer2.Vm
    int3.comp(poly_deg=3)

    # update compositions
    fer1.c[:] = int2.ci_bcc
    aus1.FDM_implicit(bc0=(1, 0, 0, int1.ci_fcc),
                      bcn=(1, 0, 0, int2.ci_fcc))
    # fer2.c[:] = np.linspace(int2.ci_bcc, int3.ci_bcc, fer2.n)
    fer2.c[:] = np.linspace(int2.ci_bcc, int3.ci_bcc, fer2.n)
    aus2.FDM_implicit(bc0=(1, 0, 0, int3.ci_fcc))

    # update position of interfaces and interpolate compositions
    fer1.update_grid(i, vn=int1.v)
    aus1.update_grid(i, v0=int1.v, vn=int2.v)
    fer2.update_grid(i, v0=int2.v, vn=int3.v)
    aus2.update_grid(i, v0=int3.v)

    log.printit(i, criteria=lambda i: (i+1) % each == 0)

log.close()

log.save_cprofiles()
log.save_properties('cavg')
log.save_properties('ci*')
log.save_properties('s*')
