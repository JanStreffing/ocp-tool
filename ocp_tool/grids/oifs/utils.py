import configparser
import re
from collections import namedtuple


def parse_griddes(griddes_string):
    """Helper function to provide grid descriptions as returned by 'cdo griddes'
    in a form that can be used to add new grids hardcoded in
    ocp_tool.grids.oifs.
    To use this function, pass it the output of
       cdo -s griddes <ICMGG????INIT>
    (as a string). The resulting dictionary contains one or more grid
    descriptions (depending on the number of grids in the ICMGG file). The grid
    description of interest can be pretty printed with, for example
        pp = pprint.PrettyPrinter(width=80, compact=True)
        pp.pprint(d)
    to be close to what is needed in the hard-coded grid description files.
    To further aid the process, a Jupyter Notebook is provided together with
    this code."""

    # This is what `cdo -s griddes` supposedly returns from ICMGG???INIT files.
    # The callable on the right hand sides are used to cast the values to their
    # correct types.
    griddes_types = dict(
        gridtype=str,
        gridsize=int,
        xsize=int,
        ysize=int,
        numlpe=int,
        xname=str,
        xlongname=lambda s: str(s[1:-1]),  # strip quotes
        xunits=lambda s: str(s[1:-1]),
        yname=str,
        ylongname=lambda s: str(s[1:-1]),
        yunits=lambda s: str(s[1:-1]),
        xvals=lambda s: [*map(float, s.split())],  # list of floats
        yvals=lambda s: [*map(float, s.split())],
        reducedpoints=lambda s: [*map(int, s.split())],  # list of ints
        rowlon=lambda s: [*map(int, s.split())],  # list of ints
    )

    cfg = configparser.ConfigParser()
    cfg.read_string(
        re.sub(
            r"#\s*\n#\s+(.*)\n#\s*\n",  # Matches the cdo griddes headers
            r"[\1]\n",  # and replaces by configparser headers
            griddes_string,
            re.MULTILINE,
        )
    )

    griddes = dict()
    for sec in cfg.sections():
        griddes[sec] = dict()
        for key, val in cfg[sec].items():
            try:
                griddes[sec][key] = griddes_types[key](val)
            except KeyError:
                griddes[sec][key] = val

    return griddes


def namedtuple_from_dict(name, dict_):
    return namedtuple(name, dict_.keys())(**dict_)
