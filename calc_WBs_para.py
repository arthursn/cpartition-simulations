if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        import numpy as np
        from cpartition import BCC, FCC, WBs, x2wp
        from scipy.optimize import bisect
        from scipy.interpolate import interp1d

        y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
                 Si=5.02504411E-2, Fe=9.4414085022e-1)

        for T in sys.argv[1:]:
            try:
                T = float(T)
            except:
                print('{} is not a valid temperature'.format(T))
                continue

            tdata_fcc = 'thermo/FoFo/TCFE8/{:.0f}-fcc.txt'.format(T)
            tdata_bcc = 'thermo/FoFo/TCFE8/{:.0f}-bcc.txt'.format(T)
            # tdata_fcc = 'thermo/FoFo/TCFE0/375-FCC.TXT'
            # tdata_bcc = 'thermo/FoFo/TCFE0/375-BCC.TXT'
            try:
                aust = FCC(T_C=T, tdata=tdata_fcc)
                ferr = BCC(T_C=T, tdata=tdata_bcc, E=WBs(T))
            except:
                print('Failed to initialize ferr and/or aust')
                continue

            # We want to find the zero of this function, which
            # is the difference between muZ(aust) and muZ(ferr)
            # both expressed as functions of muC
            # When g(muC) = 0, then muC and muZ are equal for
            # ferr and aust
            def g(muC): return aust.muC2muZ(muC) - ferr.muC2muZ(muC)

            # minimum and maximum values of mu_C
            lo = max(min(aust.chempot['MU(C)']), min(ferr.chempot['MU(C)']))
            hi = min(max(aust.chempot['MU(C)']), max(ferr.chempot['MU(C)']))

            muC = bisect(g, lo, hi, xtol=1e-3)
            cferr = ferr.mu2x['C'](muC)
            caust = aust.mu2x['C'](muC)

            print('{:g} oC, muC={:g} J/mol, caust={:g} wt.%'.format(T, muC, x2wp(caust, y=y)))
