if __name__ == '__main__':
    import sys
    from cpartition import FCC, x2wp

    if len(sys.argv) > 0:
        T_C = 375.
        c0 = 3.34414e-02
        y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
                 Si=5.02504411E-2, Fe=9.4414085022e-1)

        tdata_fcc = 'thermo/FoFo/TCFE8/375-fcc.txt'
        aust = FCC(T_C=T_C, tdata=tdata_fcc)

        for muC in sys.argv[1:]:
            try:
                muC = float(muC)
            except ValueError:
                print('Cannot parse {} to float'.format(muC))
            except:
                print('Unexpected error')
            else:
                cCCEtheta = x2wp(aust.mu2x['C'](muC), y=y)
                print('muC={:g}, ci_fcc={:g}'.format(muC, cCCEtheta))
