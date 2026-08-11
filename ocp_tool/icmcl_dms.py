"""Add the DMS surface-emission climatology to the OpenIFS ICMCL file.

WHY THIS EXISTS.  OpenIFS marine cloud condensation nuclei come from Genthon (1992)
sea-salt-from-wind alone (surfrad_ctl_mod.F90, LCCNO branch).  There is no biogenic
term, and MACv2-SP -- the only other live CCN path -- is anthropogenic plumes only.  So
natural marine biogenic sulphate is absent by construction, and the consequence is
measurable: the model's Southern Ocean droplet number is correct in the mean
(Nd ~ 73 cm-3, observed 60-120) but its SEASONAL PHASE IS INVERTED -- lowest in DJF
(68.5) and highest in JJA (76.1), because it follows the wind, where observations peak in
DJF with the biogenic bloom.  The measured SO shortwave CRE error is worst in exactly the
months the model has fewest droplets (DJF +10.86, SON +13.11 W/m2).

THE DATA ALREADY SHIPS AND NOTHING READS IT.  `climate.v020/<res>/month_dms` is a 12-month
DMS surface-emission climatology (paramId 210043, kg m-2 s-1) provided by ECMWF at every
resolution, already on the model's own reduced Gaussian grid.  It is not ingested: ICMCL
carries only albedo and LAI, and OpenIFS references DMS solely in the TM5/M7 chemistry
paths, which need NACTAERO /= 0.

WHY THIS BELONGS IN OCP-TOOL AND NOT mkclimr.sh.  emdms is an OCEAN field carrying a -99
sentinel over land -- the opposite convention to albedo and LAI, which are land fields.
Under a changed coastline (any paleo configuration) new ocean points need infill and new
land points need the sentinel.  ocp-tool already applies the land-sea mask to ICMCL in
paleo_input.py (step 14b); mkclimr.sh cannot.  Putting the field here is what makes it
paleo-safe.

THE READER'S CONTRACT, which sets the output ordering.  ece_updclie_climr.F90 discovers
the field list from the FIRST date and then requires every later date to repeat it in the
same order, aborting with "Inconsistent order of fields in ICMCL*INIT" otherwise.  Fields
are grouped BY DATE.  So the emdms message for each month must be interleaved into that
month's group, not appended in a block at the end.

Capacity was checked before writing this: NCLIGC is dimensioned 19, max_fields is 10, and
CLIMR is allocated from the discovered count -- so an eighth field needs no model change
to be READ.  That makes the file side testable on its own, before any physics is wired up.
"""
import eccodes as ecc

DMS_PARAM_ID = 210043
DMS_SHORTNAME = 'emdms'
DMS_LAND_SENTINEL = -99.0


def _messages(path):
    """Yield (dataDate, paramId, shortName, raw_bytes) for every message."""
    out = []
    with open(path, 'rb') as f:
        while True:
            gid = ecc.codes_grib_new_from_file(f)
            if gid is None:
                break
            rec = (ecc.codes_get(gid, 'dataDate'),
                   ecc.codes_get(gid, 'paramId'),
                   ecc.codes_get(gid, 'shortName'),
                   ecc.codes_get_message(gid))
            ecc.codes_release(gid)
            out.append(rec)
    return out


def add_dms_to_icmcl(icmcl_in, dms_file, icmcl_out, land_sea_mask=None,
                     verbose=True):
    """Interleave the DMS climatology into ICMCL, one message per date group.

    land_sea_mask, if given, is a 1-D array on the same reduced Gaussian grid with
    1 = land.  It is used to RE-APPLY the sentinel consistently with the mask
    ocp-tool has generated, which is the paleo-safe path: points that are ocean in
    the new mask but carry the sentinel are filled from the field's own ocean mean,
    and points that are land in the new mask are forced to the sentinel regardless
    of what the source climatology said.

    Returns the number of messages written.
    """
    icmcl = _messages(icmcl_in)
    dms = _messages(dms_file)

    dms_by_date = {}
    for date, pid, name, raw in dms:
        if pid != DMS_PARAM_ID:
            continue
        if date in dms_by_date:
            raise RuntimeError(f'duplicate DMS message for date {date}')
        dms_by_date[date] = raw
    if not dms_by_date:
        raise RuntimeError(f'no paramId {DMS_PARAM_ID} messages in {dms_file}')

    icmcl_dates = []
    for date, _, _, _ in icmcl:
        if date not in icmcl_dates:
            icmcl_dates.append(date)
    missing = [d for d in icmcl_dates if d not in dms_by_date]
    if missing:
        raise RuntimeError(
            f'ICMCL has dates with no DMS counterpart: {missing}. The reader groups '
            'by dataDate and requires an identical field list per date, so a partial '
            'merge would abort the model.')

    if any(pid == DMS_PARAM_ID for _, pid, _, _ in icmcl):
        raise RuntimeError(f'{icmcl_in} already contains paramId {DMS_PARAM_ID}')

    # write each date group followed by that date's DMS message
    n = 0
    with open(icmcl_out, 'wb') as fout:
        for i, (date, pid, name, raw) in enumerate(icmcl):
            fout.write(raw)
            n += 1
            last_of_group = (i == len(icmcl) - 1) or (icmcl[i + 1][0] != date)
            if last_of_group:
                raw_dms = dms_by_date[date]
                if land_sea_mask is not None:
                    raw_dms = _remask(raw_dms, land_sea_mask)
                fout.write(raw_dms)
                n += 1

    if verbose:
        per_date = n // len(icmcl_dates)
        print(f'  ICMCL + DMS: {len(icmcl)} -> {n} messages '
              f'({len(icmcl_dates)} dates x {per_date} fields)')
    return n


def _remask(raw, land_sea_mask):
    """Force the land sentinel to follow the supplied mask (paleo-safe)."""
    gid = ecc.codes_new_from_message(raw)
    try:
        vals = ecc.codes_get_values(gid)
        if len(vals) != len(land_sea_mask):
            raise RuntimeError(
                f'mask length {len(land_sea_mask)} does not match DMS field '
                f'length {len(vals)}')
        is_land = land_sea_mask >= 0.5
        ocean_vals = vals[(~is_land) & (vals > DMS_LAND_SENTINEL / 2)]
        fill = float(ocean_vals.mean()) if ocean_vals.size else 0.0
        # new ocean (was sentinel, now sea) -> ocean mean; any land -> sentinel
        vals = vals.copy()
        vals[(~is_land) & (vals <= DMS_LAND_SENTINEL / 2)] = fill
        vals[is_land] = DMS_LAND_SENTINEL
        ecc.codes_set_values(gid, vals)
        return ecc.codes_get_message(gid)
    finally:
        ecc.codes_release(gid)


def verify(icmcl_out, expect_fields=8, verbose=True):
    """Check the written file against the reader's contract.

    Reproduces what ece_updclie_climr.F90 (setup_toc) does: read the field order from
    the first date, then require every later date to match it exactly.  Catching a
    violation here is seconds; catching it in the model is a failed run.
    """
    msgs = _messages(icmcl_out)
    order, dates, cur, curdate = None, 0, [], None
    for date, pid, name, _ in msgs:
        if date != curdate:
            if cur:
                if order is None:
                    order = cur
                elif cur != order:
                    raise RuntimeError(
                        f'field order differs at date {curdate}: {cur} != {order} '
                        '-- the model would abort with "Inconsistent order of fields"')
                dates += 1
            cur, curdate = [], date
        cur.append(pid)
    if cur:
        if order is not None and cur != order:
            raise RuntimeError(f'field order differs at final date {curdate}')
        order = order or cur
        dates += 1

    if len(order) != expect_fields:
        raise RuntimeError(f'expected {expect_fields} fields per date, got {len(order)}')
    if DMS_PARAM_ID not in order:
        raise RuntimeError(f'paramId {DMS_PARAM_ID} absent from the written file')
    if len(order) > 10:
        raise RuntimeError('more than max_fields=10 per date; the model would abort')

    if verbose:
        print(f'  verify: {dates} dates x {len(order)} fields, order {order} -- OK')
        print(f'          within max_fields=10 and NCLIGC(19)')
    return order
