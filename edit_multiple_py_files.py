hasdjk = {'coupled_FoFo_375_CCEortho.py': 10,
          'coupled_FoFo_375_CCEpara.py': 100,
          'coupled_FoFo_375_CCE.py': 100,
          'coupled_FoFo_375_mu16e3.py': 10,
          'coupled_FoFo_375_mu20e3.py': 10,
          'coupled_FoFo_375_mu23e3.py': 10,
          'coupled_FoFo_375_mu30e3.py': 10,
          'mart_FoFo_CCEpara.py': 100,
          'mart_FoFo_CCE.py': 1,
          'mart_FoFo_CCE.py': 100,
          'mart_FoFo_mu20e3.py': 100,
          'mart_FoFo_mu23e3.py': 100,
          'mart_FoFo_mu30e3.py': 100}

from fnmatch import fnmatch

for fname, each in hasdjk.items():
    print(fname)
    f = open(fname, 'r')
    contents = f.readlines()
    f.close()

    for idx, line in enumerate(contents):
        if fnmatch(line, 't = (np.arange*'):
            if fnmatch(contents[idx+1], 'each*'):
                print(contents[idx] + contents[idx+1])
            else:
                # print(idx, contents[idx])
                contents.insert(idx+1, 'each = {}\n'.format(each))
                break

    contents = ''.join(contents)
    f = open(fname, 'w')
    f.write(contents)
    f.close()

    # break
