def lookup_option(option, args, dtype=None, vallist=[], multi=False):
    """
    Get options in list of arguments
    """
    if option in args:
        idx = args.index(option)
        args.pop(idx)
        if dtype:
            try:
                if multi:
                    while True:
                        try:
                            val = args[idx]
                            val = dtype(val)
                        except:
                            break
                        else:
                            vallist += [val]
                            args.pop(idx)
                else:
                    val = args.pop(idx)
                    val = dtype(val)
                    vallist += [val]
            except IndexError:
                print('No argument provided')
            except ValueError:
                print('Failed parsing {} as {}'.format(val, dtype))
            except:
                print('Unexpected error')
        else:
            vallist += [True]

    while option in args:
        vallist, args = lookup_option(option, args, dtype, vallist)

    return vallist, args


def filter_number(s, dtype=int):
    if s == '':
        s = None
    else:
        try:
            s = dtype(s)
        except:
            raise
    return s


def split_string(string, dtype=int, splitchar=':'):
    spt = []
    try:
        spt = list(map(lambda s: filter_number(
            s, dtype), string.split(splitchar)))
    except:
        print('Failed parsing {}'.format(string))
    return spt
