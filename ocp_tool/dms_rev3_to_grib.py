"""Convert the DMS-Rev3 climatology onto the OpenIFS reduced Gaussian grid.

WHY REPLACE THE SHIPPED FIELD.  `climate.v020/<res>/month_dms` is Lana et al. (2011).
DMS-Rev3 (Hulswar et al. 2022, ESSD 14, 2963-2987) is the third revision and rests on
865,109 observations against L11's 47,313 -- an ~18-fold increase -- with revised
filtering, dynamic biogeochemical province boundaries and a new smoothing algorithm.
Globally it is only 4 % lower, but the paper is explicit that "the largest changes are
observed in high concentration regions such as the polar oceans", and that is exactly
where this campaign is looking.

Measured on the Southern Ocean band 65-45S, L11 against Rev3:

    annual   1.93 -> 2.21   (+14 %)
    DJF      4.23 -> 4.28   ( +1 %)
    SON      1.34 -> 2.00   (+49 %)
    November 2.18 -> 4.18   (+92 %)
    global max 53.25 -> 26.56   (-50 %, the "reduced patchiness" of the paper)

So Rev3 does not raise the summer peak, it makes it BROADER -- the bloom starts a month
earlier.  The model's SO shortwave CRE error is worst in SON (-7.03 pp) and DJF (-6.94),
so the upgrade lands in the season the old climatology was weakest.

UNITS, DERIVED NOT READ.  The Rev3 netCDF carries NO attributes whatsoever, so the units
had to be established by physical constraint: integrating DMSflux globally gives 28.19
Tg S/yr (DMSfluxSI 28.36) against a published 27.1-28.1, which fixes the flux as
mol m-2 d-1 and the concentration as nmol/L.  Any other assumption is out by orders of
magnitude.  Only the CONCENTRATION is used here -- the shipped flux is computed with
Rev3's own wind climatology, whereas the model must form the flux online from its own
wind, which is quadratic in speed.  The shipped flux is kept as a validation target.

METHOD, and why nearest-neighbour.  Rev3 is 1 deg x 1 deg; TCO95's reduced Gaussian is
~0.9 deg, i.e. comparable resolution, so nearest-valid-neighbour transfers the field
without inventing structure or bleeding land NaNs into coastal ocean.  Bilinear would
need valid-weight renormalisation at every coast for no gain at matched resolution.

THE TEMPLATE IS THE SHIPPED FILE.  Values are overwritten into the existing month_dms
messages rather than building GRIB from scratch, so grid, dataDate, paramId and every
header stay byte-identical to what the reader already accepts.  The land sentinel is
taken from the template, so the model's own land mask is preserved exactly; where Rev3
has no data at an ocean point of the template (sea ice, unsampled), the L11 value is
kept and counted, because a hole would otherwise become a sentinel and be read as land.
"""
import numpy as np

DMS_PARAM_ID = 210043
DMS_LAND_SENTINEL = -99.0


def rev3_to_template_grib(rev3_nc, template_grib, out_grib, verbose=True):
    """Write a month_dms-style GRIB carrying DMS-Rev3 concentrations.

    rev3_nc        DMS_REV3.nc as distributed (Mendeley hyn62spny2)
    template_grib  the shipped climate.v020/<res>/month_dms for this resolution
    out_grib       destination

    Returns a dict of per-month diagnostics.
    """
    import eccodes as ecc
    import xarray as xr
    from scipy.spatial import cKDTree

    ds = xr.open_dataset(rev3_nc, decode_cf=False)
    rlat = ds['lat'].values                     # 89.5 .. -89.5
    rlon = ds['lon'].values                     # -179.5 .. 179.5
    conc = ds['DMS'].values.astype('float64')   # (month, lon, lat)
    conc[conc > 1e30] = np.nan
    if conc.shape[0] != 12:
        raise RuntimeError(f'expected 12 months, got {conc.shape[0]}')

    # build one KD-tree of VALID Rev3 cells per month, in 3-D so the dateline and
    # the poles cannot produce a spurious "nearest" neighbour
    LO, LA = np.meshgrid(rlon, rlat, indexing='ij')

    def xyz(lon_deg, lat_deg):
        la = np.deg2rad(lat_deg); lo = np.deg2rad(lon_deg)
        return np.stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)], -1)

    diag = {}
    fin = open(template_grib, 'rb')
    fout = open(out_grib, 'wb')
    try:
        imonth = 0
        while True:
            gid = ecc.codes_grib_new_from_file(fin)
            if gid is None:
                break
            try:
                pid = ecc.codes_get(gid, 'paramId')
                if pid != DMS_PARAM_ID:
                    raise RuntimeError(f'template message {imonth} has paramId {pid}')
                tvals = ecc.codes_get_values(gid)
                tlat = ecc.codes_get_array(gid, 'latitudes')
                tlon = ecc.codes_get_array(gid, 'longitudes')

                cm = conc[imonth]
                ok = np.isfinite(cm)
                tree = cKDTree(xyz(LO[ok], LA[ok]))
                _, idx = tree.query(xyz(np.where(tlon > 180, tlon-360, tlon), tlat))
                new = cm[ok][idx]

                is_land = tvals <= DMS_LAND_SENTINEL/2       # template's own mask
                # keep the template's land convention exactly
                new[is_land] = DMS_LAND_SENTINEL
                # a Rev3 hole at a template OCEAN point would read as land: keep L11 there
                holes = (~is_land) & (~np.isfinite(new))
                new[holes] = tvals[holes]

                ecc.codes_set_values(gid, new)
                fout.write(ecc.codes_get_message(gid))

                sea = ~is_land
                diag[imonth+1] = dict(
                    template_mean=float(tvals[sea].mean()),
                    rev3_mean=float(new[sea].mean()),
                    holes_filled=int(holes.sum()),
                    land=int(is_land.sum()))
                imonth += 1
            finally:
                ecc.codes_release(gid)
    finally:
        fin.close(); fout.close()

    if imonth != 12:
        raise RuntimeError(f'template had {imonth} messages, expected 12')
    if verbose:
        print(f'  DMS-Rev3 -> {out_grib}')
        print(f'  {"month":>5s} {"L11 mean":>9s} {"Rev3 mean":>10s} {"change":>8s} '
              f'{"holes":>6s}')
        for m in range(1, 13):
            d = diag[m]
            ch = 100*(d['rev3_mean']-d['template_mean'])/d['template_mean']
            print(f'  {m:5d} {d["template_mean"]:9.3f} {d["rev3_mean"]:10.3f} '
                  f'{ch:+7.1f}% {d["holes_filled"]:6d}')
        print(f'  land points held at the sentinel: {diag[1]["land"]}')
    return diag
