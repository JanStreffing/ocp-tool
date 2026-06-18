"""Regenerate the full atmosphere-side coupling set for a new FESOM submesh.

On a dynamic-mesh coupling leg the ocean coastline moves with the ice sheet, so
*every* land-sea-mask-derived product has to be rebuilt from the SAME new mask,
or the atmosphere and ocean disagree about the coastline:

  * the OIFS gridpoint init  ICMGG<expid>INIT  (lsm + slt fields),
  * the OASIS grids/masks/areas (A096 atmosphere, feom ocean, RnfO runoff),
  * the runoff-mapper land-sea mask,
  * the LPJ-GUESS soil-type (slt) input.

The original feom-only ``oasis_weights.regenerate_for_mesh`` only refreshed the
ocean side, leaving the atmosphere mask + ICMGG stale from chunk 2 on. This
module instead drives the *full* ocp-tool pipeline (run_ocp_tool.py) against the
submesh, so a single ``process_land_sea_mask`` produces a self-consistent set,
and then collects the products the coupling needs to stage for the next leg.

CLI:
    python -m ocp_tool.dynamic_regen \
        --mesh-dir   <submesh dir with nod2d.out/elem2d.out> \
        --grid-name  <tag, e.g. feom or submesh_19020101> \
        --template   configs/TCO95_CORE2.yaml \
        --out-dir    <couple_dir>/oasis_regen \
        [--ocp-tool-dir <repo>] [--python <interp>] [--no-rmp]

After it runs, ``--out-dir`` holds grids/masks/areas(+rmp_*).nc and the modified
ICMGG<expid>INIT_<grid-name>; the printed JSON lists the staged paths.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


def _build_config(template, mesh_dir, grid_name, ocp_tool_dir, generate_rmp):
    """Return a config dict derived from ``template`` for the submesh."""
    raw = yaml.safe_load(Path(template).read_text())

    ocean = raw.setdefault("ocean", {})
    ocean["grid_name"] = grid_name
    # read_fesom_grid_polygon builds mesh.nc from the ASCII files in the parent
    # of mesh_file when it is missing / force_overwrite_griddes is set, so point
    # mesh_file inside the submesh dir and force a rebuild for the new node set.
    ocean["mesh_file"] = str(Path(mesh_dir) / "mesh.nc")
    ocean["force_overwrite_griddes"] = True

    raw.setdefault("options", {})["generate_rmp"] = bool(generate_rmp)
    # Pin root_dir so the input/ dirs resolve regardless of where we write the
    # generated config.
    raw.setdefault("paths", {})["root_dir"] = str(Path(ocp_tool_dir).resolve())
    return raw


def regenerate_for_submesh(
    mesh_dir,
    grid_name,
    out_dir,
    template,
    ocp_tool_dir=None,
    python_exe=None,
    generate_rmp=True,
):
    """Run the full ocp-tool pipeline for ``mesh_dir`` and stage its products.

    Returns a dict with the staged ``oasis_dir`` and ``icmgg`` paths.
    """
    mesh_dir = Path(mesh_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ocp_tool_dir = Path(ocp_tool_dir) if ocp_tool_dir else Path(__file__).resolve().parent.parent
    python_exe = python_exe or sys.executable

    raw = _build_config(template, mesh_dir, grid_name, ocp_tool_dir, generate_rmp)
    cfg_path = out_dir / f"ocp_config_{grid_name}.yaml"
    cfg_path.write_text(yaml.safe_dump(raw, sort_keys=False))

    # Run the pipeline. run_ocp_tool.py lives at the repo root and imports the
    # ocp_tool package, so launch it with cwd at the repo root.
    cmd = [python_exe, "run_ocp_tool.py", str(cfg_path)]
    print(f"[dynamic_regen] {' '.join(cmd)} (cwd={ocp_tool_dir})", flush=True)
    subprocess.run(cmd, cwd=str(ocp_tool_dir), check=True)

    # Collect products from output/TCO{res}_{grid_name}/.
    resolution = raw["atmosphere"]["resolution_list"][0]
    trunc = raw["atmosphere"]["truncation_type"]
    prefix = "TCO" if trunc == "cubic-octahedral" else "TL"
    out_base = ocp_tool_dir / "output" / f"{prefix}{resolution}_{grid_name}"
    oasis_src = out_base / "oasis_mct3_input"
    icmgg_name = f"ICMGG{raw['atmosphere']['experiment_name']}INIT_{grid_name}"
    icmgg_src = out_base / "openifs_input_modified" / icmgg_name

    staged = {"oasis_files": [], "icmgg": None, "runoff_files": [],
              "slt_files": [], "config": str(cfg_path)}

    # OASIS grids/masks/areas (+ rmp_* when generate_rmp) + CO2/veg restarts.
    for f in sorted(oasis_src.glob("*.nc")):
        dst = out_dir / f.name
        shutil.copy2(f, dst)
        staged["oasis_files"].append(str(dst))

    # Modified OIFS gridpoint init (new lsm/slt).
    if icmgg_src.exists():
        dst = out_dir / icmgg_name
        shutil.copy2(icmgg_src, dst)
        staged["icmgg"] = str(dst)
    else:
        print(f"[dynamic_regen] WARNING: modified ICMGG not found at {icmgg_src}",
              file=sys.stderr)

    # Runoff-mapper land-sea mask (changes with the coastline) and LPJ-GUESS
    # soil-type input -- the other two products of the same land-sea mask.
    for f in sorted((out_base / "runoff_map_modified").glob("*.nc")):
        dst = out_dir / f.name
        shutil.copy2(f, dst)
        staged["runoff_files"].append(str(dst))
    for f in sorted((out_base / "lpj-guess").glob("*.nc")):
        dst = out_dir / f.name
        shutil.copy2(f, dst)
        staged["slt_files"].append(str(dst))

    print(json.dumps(staged, indent=2))
    return staged


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mesh-dir", required=True,
                   help="submesh dir holding nod2d.out/elem2d.out")
    p.add_argument("--grid-name", required=True,
                   help="ocean grid tag (output dir + ICMGG suffix)")
    p.add_argument("--out-dir", required=True,
                   help="where to stage grids/masks/areas/rmp + modified ICMGG")
    p.add_argument("--template", required=True,
                   help="template ocp-tool config (e.g. configs/TCO95_CORE2.yaml)")
    p.add_argument("--ocp-tool-dir", default=None,
                   help="ocp-tool repo root (default: inferred from this file)")
    p.add_argument("--python", default=None,
                   help="interpreter to run run_ocp_tool.py (default: this one)")
    p.add_argument("--no-rmp", action="store_true",
                   help="skip OASIS rmp weight generation (LSM/ICMGG/grids only)")
    args = p.parse_args(argv)

    regenerate_for_submesh(
        mesh_dir=args.mesh_dir,
        grid_name=args.grid_name,
        out_dir=args.out_dir,
        template=args.template,
        ocp_tool_dir=args.ocp_tool_dir,
        python_exe=args.python,
        generate_rmp=not args.no_rmp,
    )


if __name__ == "__main__":
    main()
