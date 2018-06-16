#!/usr/bin/env python

import numpy as np
import time

import os
import sys
from cpartition import *

basename = os.path.basename(__file__).replace('.py', '')

# wc0 = 0.76e-2
c0 = 3.34414e-02
T_C = 375.

# dt = total_time/n_time
dt = 5e-4
total_time = 10
n_time = int(total_time/dt)
dt = total_time/n_time
t = (np.arange(n_time) + 1)*dt

tdata_fcc = 'thermo/FoFo/TCFE8/375-fcc.txt'
tdata_bcc = 'thermo/FoFo/TCFE8/375-bcc.txt'
# tdata_fcc = 'thermo/FoFo/TCFE0/375-FCC.TXT'
# tdata_bcc = 'thermo/FoFo/TCFE0/375-BCC.TXT'

mart = BCC(T_C=T_C, dt=dt, z=np.linspace(-1.16, -.66, 50), c0=c0,
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

fer1.c[:] = int3.CCE(c0)
fer2.c[:] = int4.CCE(c0)

j, fer1_diss = -1, False

log = SimulationLog(basename)
log.set_domains([('mart', mart), ('aus1', aus1),
                 ('fer1', fer1), ('aus2', aus2), ('fer2', fer2)])
log.set_interfaces([('int1', int1), ('int2', int2),
                    ('int4', int4), ('int4', int4)])
log.set_conditions(c0, T_C, total_time, n_time)
log.initialize(each=100, flush=False)

for i in range(n_time):
    try:

        if not fer1_diss:
            # interface velocities at the mobile interfaces
            int1.comp(poly_deg=2)
            int2.v = 1e6*int2.chem_driving_force()*int2.M()/fer1.Vm
            int2.comp(poly_deg=2)
            int3.v = 1e6*int3.chem_driving_force()*int3.M()/fer1.Vm
            int3.comp(poly_deg=2)
            int4.v = 1e6*int4.chem_driving_force()*int4.M()/fer2.Vm
            int4.comp(poly_deg=2)

            pos0 = fer1.z[0] + int2.v*dt
            posn = fer1.z[-1] + int3.v*dt

            if pos0 < posn:
                mart.FDM_implicit(bcn=(1., 0, 0, int1.ci_bcc))
                aus1.FDM_implicit(bc0=(1, 0, 0, int1.ci_fcc),
                                  bcn=(1, 0, 0, int2.ci_fcc))
                fer1.c[:] = np.linspace(int2.ci_bcc, int3.ci_bcc, fer1.n)
                aus2.FDM_implicit(bc0=(1, 0, 0, int3.ci_fcc),
                                  bcn=(1, 0, 0, int4.ci_fcc))
                fer2.c.fill(int4.ci_bcc)

                mart.update_grid()
                aus1.update_grid(vn=int2.v)
                fer1.update_grid(v0=int2.v, vn=int3.v)
                aus2.update_grid(v0=int3.v, vn=int4.v)
                fer2.update_grid(v0=int4.v)
            else:
                # Initialize new configuration
                fer1_diss = True
                z, c, cavg = merge_domains(aus1, aus2)
                dz = np.min([aus1.dz, aus2.dz])
                n = int(np.abs((z[-1] - z[0])/dz))
                print(i+1, n)

                aus1.z = np.linspace(z[0], z[-1], n)
                aus1.c = interp1d(z, c)(aus1.z)
                aus1.initialize_grid(reset=False)

                int1 = Interface(domain1=mart, domain2=aus1, type_int='fixed.fluxes')
                int2 = Interface(domain1=aus1, domain2=fer2, type_int='mobile.mmode')

                fer1.deactivate()
                aus2.deactivate()

        if fer1_diss:
            int1.comp(poly_deg=2)
            int2.v = 1e6*int2.chem_driving_force()*int2.M()/fer2.Vm
            int2.comp(poly_deg=2)

            mart.FDM_implicit(bcn=(1., 0, 0, int1.ci_bcc))
            aus1.FDM_implicit(bc0=(1, 0, 0, int1.ci_fcc),
                              bcn=(1, 0, 0, int2.ci_fcc))
            fer2.c.fill(int2.ci_bcc)

            mart.update_grid()
            aus1.update_grid(vn=int2.v)
            aus2.update_grid()
            fer1.update_grid()
            fer2.update_grid(v0=int2.v)

            j += 1

    except:
        print(i+1, j)
        raise

    log.print(i)

log.close()

log.save_cprofiles()
log.save_properties('cavg')
log.save_properties('ci*')
log.save_properties('s*')
