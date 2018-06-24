import numpy as np
import time

import os
import sys
from cpartition import *

basename = os.path.basename(__file__).replace('.py', '')

# wc0 = 0.76e-2
c0 = 3.34414e-02
T_C = 375.

control_itsteps = ControlIterationSteps([5e-4, 5e-3, 5e-2], [0, 1, 150, 1000])
total_time = control_itsteps.total_time
n_time = control_itsteps.ntime
dt = control_itsteps.dt
each = 20
control_itsteps.print_summary()

tdata_fcc = 'thermo/FoFo/TCFE8/375-fcc.txt'
tdata_bcc = 'thermo/FoFo/TCFE8/375-bcc.txt'

mart = BCC(T_C=T_C, dt=dt, z=np.linspace(-1.16, -.66, 20), c0=c0,
           n_time=n_time, tdata=tdata_bcc,
           type_D='carbides', cmax_bcc=5.4e-4, c_carbide=.25)
aus1 = FCC(T_C=T_C, dt=dt, z=np.linspace(-.66, -.33, 100), c0=c0,
           n_time=n_time, tdata=tdata_fcc)
fer1 = BCC(T_C=T_C, dt=dt, z=np.linspace(-.33, -.33, 10), c0=0.,
           n_time=n_time, tdata=tdata_bcc, E=WBs(T_C))
aus2 = FCC(T_C=T_C, dt=dt, z=np.linspace(-.33, 0, 100), c0=c0,
           n_time=n_time, tdata=tdata_fcc)
fer2 = BCC(T_C=T_C, dt=dt, z=np.linspace(0, 0, 10), c0=0.,
           n_time=n_time, tdata=tdata_bcc, E=WBs(T_C))

int1 = Interface(domain1=mart, domain2=aus1, type_int='fixed.fluxes')
int2 = Interface(domain1=aus1, domain2=fer1, type_int='mobile.mmode')
int3 = Interface(domain1=fer1, domain2=aus2, type_int='mobile.mmode')
int4 = Interface(domain1=aus2, domain2=fer2, type_int='mobile.mmode')

pos_mart = mart.z[-1]
j, aus1_diss = -1, False

# fixed composition set by CCEtheta in austenite at the interface
muC = 8188.68
int1.ci_fcc = aus1.mu2x['C'](muC)
print('muC={:}, ci_fcc={:}\n'.format(muC, int1.ci_fcc))

fer1.c[:] = int3.CCE(c0)
fer2.c[:] = int4.CCE(c0)

log = SimulationLog(basename)
log.set_domains([('mart', mart), ('aus1', aus1),
                 ('fer1', fer1), ('aus2', aus2), ('fer2', fer2)])
log.set_interfaces([('int1', int1), ('int2', int2),
                    ('int4', int4), ('int4', int4)])
log.set_conditions(c0, T_C, total_time, n_time)
log.initialize(False)

for i in control_itsteps.itlist:
    if i in control_itsteps.itstepi and i > 0:
        control_itsteps.next_itstep()
        dt = control_itsteps.dt

        mart.dt = dt
        aus1.dt = dt
        fer1.dt = dt
        aus2.dt = dt
        fer2.dt = dt

    try:
        if not aus1_diss:
            # interface velocities at the mobile interfaces
            int2.v = 1e6*int2.chem_driving_force()*int2.M()/fer1.Vm
            int2.comp(poly_deg=3)
            int3.v = 1e6*int3.chem_driving_force()*int3.M()/fer1.Vm
            int3.comp(poly_deg=3)
            int4.v = 1e6*int4.chem_driving_force()*int4.M()/fer2.Vm
            int4.comp(poly_deg=3)

            pos_fer1 = fer1.z[0] + int2.v*dt
            if pos_fer1 > pos_mart:
                if aus1.r.max() > 1. and len(aus1.z) > 3:
                    n = int(len(aus1.z)/2)
                    z, c = aus1.z, aus1.c
                    aus1.z = np.linspace(z[0], z[-1], n)
                    aus1.c = interp1d(z, c)(aus1.z)
                    aus1.initialize_grid()
                    print('aus1', i+1, n)

                J, = int1.flux('fcc')
                mart.FDM_implicit(bcn=(1.5, -2., .5, -J*mart.dz/mart.D(C=mart.c[-1])))
                aus1.FDM_implicit(bc0=(1, 0, 0, int1.ci_fcc),
                                  bcn=(1, 0, 0, int2.ci_fcc))
                fer1.c[:] = np.linspace(int2.ci_bcc, int3.ci_bcc, fer1.n)
                aus2.FDM_implicit(bc0=(1, 0, 0, int3.ci_fcc),
                                  bcn=(1, 0, 0, int4.ci_fcc))
                fer2.c.fill(int4.ci_bcc)

                # Update position of interfaces and interpolate compositions
                mart.update_grid(i)
                aus1.update_grid(i, vn=int2.v)
                fer1.update_grid(i, v0=int2.v, vn=int3.v)
                aus2.update_grid(i, v0=int3.v, vn=int4.v)
                fer2.update_grid(i, v0=int4.v)
            else:
                aus1_diss = True
                # z, c, cavg = merge_domains(mart, fer1)
                # n = int(np.abs((z[-1] - z[0])/mart.dz))
                # print('fer1', i+1, n)

                # fer1.z = np.linspace(z[0], z[-1], n)
                # fer1.c = interp1d(z, c)(fer1.z)
                # fer1.initialize_grid()

                # mart.deactivate()
                aus1.deactivate()

        if aus1_diss:
            int3.v = 1e6*int3.chem_driving_force()*int3.M()/fer1.Vm
            int3.comp(poly_deg=3)
            int4.v = 1e6*int4.chem_driving_force()*int4.M()/fer2.Vm
            int4.comp(poly_deg=3)

            mart.FDM_implicit()
            # fer1.FDM_implicit(bcn=(1, 0, 0, int3.ci_bcc))
            fer1.c.fill(int3.ci_bcc)
            aus2.FDM_implicit(bc0=(1, 0, 0, int3.ci_fcc),
                              bcn=(1, 0, 0, int4.ci_fcc))
            fer2.c.fill(int4.ci_bcc)

            mart.update_grid(i)
            aus1.update_grid(i)
            fer1.update_grid(i, vn=int3.v)
            aus2.update_grid(i, v0=int3.v, vn=int4.v)
            fer2.update_grid(i, v0=int4.v)

            j += 1

    except:
        print('error', i+1, j)
        raise

    log.printit(i, criteria=lambda i: (i+1) % each == 0)

log.close()

log.save_cprofiles()
log.save_properties('cavg')
log.save_properties('ci*')
log.save_properties('s*')
