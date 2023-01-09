"""GRIB file handling with ecCodes
"""
import eccodes as ecc


def read(file, shortnames):
    """Reads all messages in a grib file, checks the 'shortName' grib key
    against the given shortname iterable and returns a dict with values for
    all matching messages. Returns None for non-existing shortnames.
    Raises RuntimeError if the same name appears in more than one grib
    message.
    """
    data = {name: None for name in shortnames}
    with open(file, 'rb') as f:
        while True:
            gid = ecc.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                name = ecc.codes_get(gid, 'shortName')
            except ecc.KeyValueNotFoundError:
                pass
            else:
                if name in shortnames:
                    if data[name] is not None:
                        raise RuntimeError(
                            f'shortName {name} found in more than '
                            'one grib messages'
                        )
                    data[name] = ecc.codes_get_values(gid)
            ecc.codes_release(gid)
    return data


def copy_modify(infile, outfile, data=None):
    """Reads a GRIB file (infile) and copies all messages to another file
    (outfile). If the 'data' dict is given, the corresponding data with
    matching shortName in the inflile is replaced.
    """
    with open(infile, 'rb') as fin, \
         open(outfile, 'wb') as fout:

        while True:
            gid = ecc.codes_grib_new_from_file(fin)
            if gid is None:
                break
            try:
                name = ecc.codes_get(gid, 'shortName')
            except ecc.KeyValueNotFoundError:
                pass
            else:
                if data is not None and name in data:
                    ecc.codes_set_values(gid, data[name])
            ecc.codes_write(gid, fout)
            ecc.codes_release(gid)
