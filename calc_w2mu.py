if __name__ == '__main__':
    import sys

    if len(sys.argv) > 0:
        from cpartition import FCC, w2x

        T_C = 375.
        c0 = 3.34414e-02
        y = dict(Cu=3.55354266E-3, Mn=2.05516602E-3,
                 Si=5.02504411E-2, Fe=9.4414085022e-1)

        tdata_fcc = 'thermo/FoFo/TCFE8/375-fcc.txt'
        aust = FCC(T_C=T_C, tdata=tdata_fcc)

        for wC in sys.argv[1:]:
            try:
                wC = float(wC)
            except ValueError:
                print('Cannot parse {} to float'.format(wC))
            except:
                print('Unexpected error')
            else:
                xC = w2x(1e-2*wC, y=y)
                muC = float(aust.x2mu['C'](xC))
                print('muC={:g}, ci_fcc={:g}'.format(muC, wC))
