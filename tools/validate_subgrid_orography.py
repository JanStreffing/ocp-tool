#!/usr/bin/env python3
"""Check the Lott and Miller subgrid orography against what ECMWF ships.

Nothing that writes an ICMGG should be trusted until this reproduces
``climate.v020/<res>/{stdgwd,isogwd,slogwd}``. Run it as

    python tools/validate_subgrid_orography.py \
        --climate-dir /work/ab0246/a270092/input/oifs-48r1/climate.v020/95_4 \
        --fine /work/ab0246/a270092/input/oifs-48r1/etopo1/ETOPO1_Bed_c_gmt4.grd

A global pass over ETOPO1 takes about 15 seconds.

With the Bed version of ETOPO1 the comparison must exclude the ice sheets,
which is what ``--exclude-ice-sheets`` does: Bed carries bedrock under
Greenland and Antarctica, and comparing that against an ECMWF field built from
the ice surface says nothing about the formulation. With the Ice version, drop
the flag.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ocp_tool.subgrid_orography import (  # noqa: E402
    FineOrography,
    compute_sso,
    format_report,
    read_target_grid,
    validate,
)


def grib_values(path, want_coords=False):
    import eccodes as ec

    with open(path, "rb") as f:
        gid = ec.codes_grib_new_from_file(f)
        try:
            vals = np.asarray(ec.codes_get_values(gid), dtype=float)
            if not want_coords:
                return vals
            lat = np.asarray(ec.codes_get_array(gid, "latitudes"), dtype=float)
            lon = np.asarray(ec.codes_get_array(gid, "longitudes"), dtype=float)
            return vals, lat, lon
        finally:
            ec.codes_release(gid)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--climate-dir", required=True, type=Path,
                   help="ECMWF climate set for this resolution, e.g. .../climate.v020/95_4")
    p.add_argument("--fine", required=True, type=Path,
                   help="fine orography netCDF (ETOPO1 gmt4 .grd)")
    p.add_argument("--variable", default="z")
    p.add_argument("--band-rows", type=int, default=600)
    p.add_argument("--exclude-ice-sheets", action="store_true",
                   help="drop Antarctica and Greenland, required with ETOPO1 Bed")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()

    target = args.climate_dir / "stdgwd"
    grid = read_target_grid(target)
    print(f"target grid: {grid.nrows} rows, {grid.size} points, "
          f"row lengths {grid.pl.min()}-{grid.pl.max()}")

    fine = FineOrography(path=args.fine, variable=args.variable,
                         band_rows=args.band_rows)
    print(f"fine orography: {args.fine.name}")
    out = compute_sso(fine, grid, verbose=args.verbose)

    lsm, lat, lon = grib_values(args.climate_dir / "lsmoro", want_coords=True)
    mask = lsm > 0.5
    if args.exclude_ice_sheets:
        lon360 = (lon + 360.0) % 360.0
        antarctica = lat < -60.0
        greenland = (lat > 59.0) & (lat < 84.0) & (lon360 > 285.0) & (lon360 < 350.0)
        mask &= ~(antarctica | greenland)
        print("excluded Antarctica and Greenland (ETOPO1 Bed is bedrock there)")

    # The box-mean orography is the control: it has nothing to do with the SSO
    # formulation, so if it does not match, the grid mapping is wrong and every
    # other number below is meaningless.
    oro = grib_values(args.climate_dir / "orog")
    d = out["mean"][mask] - oro[mask]
    corr = np.corrcoef(out["mean"][mask], oro[mask])[0, 1]
    print(f"\ncontrol, box-mean orography: mine {out['mean'][mask].mean():.1f} m, "
          f"ECMWF {oro[mask].mean():.1f} m, rms {np.sqrt((d ** 2).mean()):.1f} m, "
          f"corr {corr:.4f}")
    if corr < 0.95:
        print("  WARNING: the mean orography does not match, so the grid mapping "
              "is wrong. Fix that before reading anything below.")

    print()
    print(format_report(validate(out, args.climate_dir, land_mask=mask)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
