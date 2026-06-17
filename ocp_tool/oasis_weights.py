"""Generate OASIS3-MCT remapping weight files (rmp_*.nc) with pyOASIS.

This replaces the earlier non-working ``oasis_remap.py`` stub. The design is
ported from rdy2cpl (https://github.com/uwefladrich/rdy2cpl, U. Fladrich): run
the OASIS coupler itself, one MPI rank per *unique* coupling link, let each rank
define its source/target grids and a coupled variable pair, and call
``enddef()`` so OASIS/SCRIP computes and writes that link's weight file.

Differences from rdy2cpl:

* pyOASIS is taken from the *coupled model's own* OASIS build (the one esm-tools
  compiles with ``OASIS_WITH_PYOASIS=ON``), not a separate OASIS install. See
  ``import_pyoasis``. This avoids maintaining a second OASIS just for weights.
* Grid geometry (centres, corners, mask, area) is read back from the OASIS grid
  description files (grids.nc / masks.nc / areas.nc) that ocp-tool already
  writes for every grid -- including the ``feom`` grid produced for an arbitrary
  (e.g. dynamically remeshed) FESOM mesh by ``pyfesom2.write_fesom_oasis_files``.
  So no grid geometry is recomputed here; we only build the OASIS weights.

The module is both an importable driver (``generate_weights``) and an MPI worker
(``python -m ocp_tool.oasis_weights <run_dir>``); the driver launches the worker
under ``srun``/``mpirun`` with one task per link.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import numpy as np
from netCDF4 import Dataset


# ---------------------------------------------------------------------------
# pyOASIS loader -- use the coupled model's OASIS build
# ---------------------------------------------------------------------------
def import_pyoasis(oasis_build_path: Optional[str] = None):
    """Import and return the ``pyoasis`` module from the model's OASIS build.

    Resolution order for the OASIS build directory:
      1. explicit ``oasis_build_path`` argument,
      2. ``$OASIS_BUILD_PATH`` (the directory containing ``python/pyoasis`` and
         ``lib/liboasis.cbind.so``),
      3. ``$OASIS_DIR``.

    pyOASIS dlopen()s ``liboasis.cbind.so`` at import, so that library's
    directory must be on ``LD_LIBRARY_PATH`` before the Python process starts
    (set by the caller's environment / the esm-tools coupling env).
    """
    build = oasis_build_path or os.environ.get("OASIS_BUILD_PATH") or os.environ.get("OASIS_DIR")
    if not build:
        raise ImportError(
            "Cannot locate the OASIS build: set OASIS_BUILD_PATH (or OASIS_DIR) "
            "to the model's oasis directory (the one containing python/pyoasis "
            "and lib/liboasis.cbind.so)."
        )
    build = str(Path(build).resolve())
    pydir = os.path.join(build, "python")
    if pydir not in sys.path:
        sys.path.insert(0, pydir)

    libdir = os.path.join(build, "lib")
    if "LD_LIBRARY_PATH" not in os.environ or libdir not in os.environ["LD_LIBRARY_PATH"].split(":"):
        # Informative only: dlopen uses the loader's LD_LIBRARY_PATH, which a
        # running process cannot change for itself. Warn so the launcher fixes it.
        sys.stderr.write(
            f"[oasis_weights] warning: {libdir} not on LD_LIBRARY_PATH; "
            "pyOASIS may fail to load liboasis.cbind.so\n"
        )
    try:
        import pyoasis  # noqa: F401
    except Exception as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            f"Failed to import pyoasis from {pydir} (build={build}): {exc}"
        ) from exc
    return pyoasis


# ---------------------------------------------------------------------------
# Coupling-link description
# ---------------------------------------------------------------------------
@dataclass
class Link:
    """One OASIS coupling link == one remapping weight file rmp_<src>_to_<tgt>_<MAP>.nc."""

    source: str            # source grid name in grids.nc (e.g. "A096")
    target: str            # target grid name (e.g. "feom")
    scripr: List[str]      # SCRIPR option lines, e.g. ["GAUSWGT D SCALAR LATITUDE 1 25 0.1", "GLBPOS opt"]
    loctrans: Optional[str] = None  # optional LOCTRANS op (e.g. "AVERAGE") -- not needed for weights

    @property
    def map_name(self) -> str:
        # scripr == ["SCRIPR", "<MAP> ...", "<opt>", ...]; the map kind is the
        # first token of the line after the SCRIPR keyword.
        for line in self.scripr:
            tok = line.split()[0]
            if tok != "SCRIPR":
                return tok
        return "SCRIPR"  # GAUSWGT / BICUBIC / CONSERV / DISTWGT

    def as_dict(self):
        return {"source": self.source, "target": self.target,
                "scripr": self.scripr, "loctrans": self.loctrans}

    @staticmethod
    def from_dict(d):
        return Link(d["source"], d["target"], d["scripr"], d.get("loctrans"))


# Default awiesm3 feom links -- these reproduce the remapping methods of the
# existing TCO95/CORE2 namcouple (non-conservative, centre based). When the
# feom grid is regenerated with corners, set ``method="conserv"`` in
# generate_weights to switch the flux links to conservative remapping.
def awiesm3_feom_links(atm_grid: str = "A096", method: str = "existing") -> List[Link]:
    if method not in ("existing", "conserv"):
        raise ValueError(f"unknown method profile: {method}")
    gauswgt = ["SCRIPR", "GAUSWGT D SCALAR LATITUDE 1 25 0.1", "GLBPOS opt"]
    gauswgt_u = ["SCRIPR", "GAUSWGT U SCALAR LATITUDE 1 25 0.1"]
    bicubic = ["SCRIPR", "BICUBIC D SCALAR LATITUDE 15"]
    conserv = ["SCRIPR", "CONSERV LR SCALAR LATITUDE 1 25 0.1", "GLBPOS opt"]
    rnf_gauswgt = ["SCRIPR", "GAUSWGT LR SCALAR LATITUDE 1 25 0.1", "GLBPOS opt"]

    flux_atm_to_oce = conserv if method == "conserv" else gauswgt
    rnf_to_oce = conserv if method == "conserv" else rnf_gauswgt

    return [
        Link(atm_grid, "feom", flux_atm_to_oce),   # A096 -> feom (heat/freshwater fluxes)
        Link(atm_grid, "feom", bicubic),           # A096 -> feom (state, bicubic)
        Link("feom", atm_grid, gauswgt_u),         # feom -> A096 (SST, ice state)
        Link("RnfO", "feom", rnf_to_oce),          # runoff -> feom
    ]


# ---------------------------------------------------------------------------
# Grid geometry read back from the OASIS description files
# ---------------------------------------------------------------------------
def _read_var(path: Path, name: str):
    if not path.exists():
        return None
    with Dataset(path) as nc:
        if name not in nc.variables:
            return None
        return np.array(nc.variables[name][:])


class OasisFileGrid:
    """A grid whose geometry is read back from grids.nc / masks.nc / areas.nc.

    OASIS stores ``<name>.lon(y, x)`` etc.; pyOASIS expects (nx, ny) Fortran
    order, so arrays are transposed on read. Corners (``<name>.clo/.cla``,
    shape (crn, y, x)) are optional -- present for grids that support
    conservative remapping (e.g. feom written with corners by pyfesom2).
    """

    def __init__(self, name: str, oasis_dir):
        self.name = name
        d = Path(oasis_dir)
        lon = _read_var(d / "grids.nc", f"{name}.lon")
        lat = _read_var(d / "grids.nc", f"{name}.lat")
        if lon is None or lat is None:
            raise KeyError(f"grid '{name}' not found in {d/'grids.nc'}")
        # (y, x) -> (x, y)
        self.center_longitudes = np.ascontiguousarray(lon.T, dtype="float64")
        self.center_latitudes = np.ascontiguousarray(lat.T, dtype="float64")
        self.shape = self.center_longitudes.shape  # (nx, ny)
        self.size = int(self.center_longitudes.size)

        clo = _read_var(d / "grids.nc", f"{name}.clo")
        cla = _read_var(d / "grids.nc", f"{name}.cla")
        if clo is not None and cla is not None:
            # (crn, y, x) -> (x, y, crn)
            self.corner_longitudes = np.ascontiguousarray(np.transpose(clo, (2, 1, 0)), dtype="float64")
            self.corner_latitudes = np.ascontiguousarray(np.transpose(cla, (2, 1, 0)), dtype="float64")
        else:
            self.corner_longitudes = None
            self.corner_latitudes = None

        srf = _read_var(d / "areas.nc", f"{name}.srf")
        self.areas = np.ascontiguousarray(srf.T, dtype="float64") if srf is not None else None

        msk = _read_var(d / "masks.nc", f"{name}.msk")
        # OASIS convention: 1 = masked (excluded), 0 = active. Default all active.
        if msk is not None:
            self.mask = np.ascontiguousarray(msk.T, dtype="int32")
        else:
            self.mask = np.zeros(self.shape, dtype="int32")

    def write_pyoasis(self, pyoasis):
        """Declare this grid to OASIS (writes into grids.nc/masks.nc/areas.nc)."""
        partition = pyoasis.SerialPartition(self.size)
        grid = pyoasis.Grid(self.name, self.shape[0], self.shape[1],
                            self.center_longitudes, self.center_latitudes,
                            partition)
        if self.corner_longitudes is not None:
            grid.set_corners(self.corner_longitudes, self.corner_latitudes)
        if self.areas is not None:
            grid.set_area(self.areas)
        grid.set_mask(self.mask)
        grid.write()
        return partition


# ---------------------------------------------------------------------------
# namcouple writer (minimal -- only what OASIS needs to compute the weights)
# ---------------------------------------------------------------------------
def write_namcouple(run_dir: Path, links: List[Link], runtime: int = 3600):
    lines = [
        "# namcouple generated by ocp_tool.oasis_weights for weight generation only",
        "$NFIELDS", f"  {len(links)}", "$END",
        "$RUNTIME", f"  {runtime}", "$END",
        "$NLOGPRT", "  1 0", "$END",
        "$NNOREST", "  False", "$END",
        "$STRINGS", "",
    ]
    for i, lk in enumerate(links):
        src_grid = OasisFileGrid(lk.source, run_dir)
        tgt_grid = OasisFileGrid(lk.target, run_dir)
        trans = ([lk.loctrans] if lk.loctrans else []) + lk.scripr
        ntrans = (1 if lk.loctrans else 0) + 1  # LOCTRANS (optional) + SCRIPR
        lines.append(f"# link {i}: {lk.source} -> {lk.target} ({lk.map_name})")
        lines.append(f" VAR_{i:02d}_S VAR_{i:02d}_T 1 {runtime} {ntrans} none EXPORTED")
        lines.append(f" {lk.source} {lk.target}")
        lines.append(" P 0 P 0")
        if lk.loctrans:
            lines.append(f" {lk.loctrans} {lk.scripr[0]}")
            for opt in lk.scripr[1:]:
                lines.append(f" {opt}")
        else:
            lines.append(f" {lk.scripr[0]}")
            for opt in lk.scripr[1:]:
                lines.append(f" {opt}")
        lines.append("")
    lines.append("$END")
    (run_dir / "namcouple").write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# MPI worker -- one rank per link, runs OASIS enddef to write the weights
# ---------------------------------------------------------------------------
def run_worker(run_dir):
    """Executed under MPI (one rank per link). Builds grids + var pair, enddef."""
    from mpi4py import MPI

    run_dir = Path(run_dir)
    pyoasis = import_pyoasis()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    links = [Link.from_dict(d) for d in json.loads((run_dir / "links.json").read_text())]
    if size < len(links):
        raise RuntimeError(f"need >= {len(links)} MPI ranks, got {size}")

    os.chdir(run_dir)

    if rank < len(links):
        lk = links[rank]
        comp = pyoasis.Component(f"ocpw{rank:02d}")
        # Read original geometry (done by all coupled ranks) before any grid
        # write truncates grids.nc; barrier guarantees the ordering.
        src = OasisFileGrid(lk.source, run_dir)
        tgt = OasisFileGrid(lk.target, run_dir)
        comm.Barrier()
        src_part = src.write_pyoasis(pyoasis)
        tgt_part = tgt.write_pyoasis(pyoasis)
        pyoasis.Var(f"VAR_{rank:02d}_S", src_part, pyoasis.OASIS.OUT)
        pyoasis.Var(f"VAR_{rank:02d}_T", tgt_part, pyoasis.OASIS.IN)
        comp.enddef()
        del comp
    else:
        # surplus ranks: join the coupling so the component MPI split is balanced
        comm.Barrier()
        pyoasis.Component(f"ocpx{rank:02d}", coupled=False)


# ---------------------------------------------------------------------------
# Driver -- prepare run dir + namcouple, launch the MPI worker
# ---------------------------------------------------------------------------
def generate_weights(
    oasis_dir,
    links: Optional[List[Link]] = None,
    *,
    atm_grid: str = "A096",
    method: str = "existing",
    launcher: Optional[List[str]] = None,
    oasis_build_path: Optional[str] = None,
):
    """Generate rmp_*.nc weight files in ``oasis_dir``.

    ``oasis_dir`` must already contain grids.nc / masks.nc / areas.nc with all
    grids referenced by ``links`` (ocp-tool writes these; the feom grid for the
    current FESOM mesh must be present).

    ``links`` defaults to the awiesm3 feom links (see ``awiesm3_feom_links``).
    ``method`` selects the remapping profile ("existing" or "conserv").
    ``launcher`` overrides the MPI launch command (default: srun / mpirun with
    one task per link).
    """
    oasis_dir = Path(oasis_dir).resolve()
    if links is None:
        links = awiesm3_feom_links(atm_grid=atm_grid, method=method)

    write_namcouple(oasis_dir, links)
    (oasis_dir / "links.json").write_text(json.dumps([lk.as_dict() for lk in links]))

    nlinks = len(links)
    env = dict(os.environ)
    build = oasis_build_path or env.get("OASIS_BUILD_PATH") or env.get("OASIS_DIR")
    if build:
        env["OASIS_BUILD_PATH"] = str(Path(build).resolve())
        libdir = str(Path(build).resolve() / "lib")
        env["LD_LIBRARY_PATH"] = libdir + ":" + env.get("LD_LIBRARY_PATH", "")

    if launcher is None:
        if subprocess.call(["bash", "-lc", "command -v srun >/dev/null 2>&1"]) == 0:
            launcher = ["srun", "-n", str(nlinks), "--cpu-bind=none"]
        else:
            launcher = ["mpirun", "-n", str(nlinks)]

    cmd = launcher + [sys.executable, "-m", "ocp_tool.oasis_weights", str(oasis_dir)]
    print(f"[oasis_weights] launching: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=oasis_dir, env=env)

    produced = sorted(p.name for p in oasis_dir.glob("rmp_*.nc"))
    print(f"[oasis_weights] produced {len(produced)} weight files: {produced}")
    return [oasis_dir / n for n in produced]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("usage: python -m ocp_tool.oasis_weights <run_dir>\n")
        sys.exit(2)
    run_worker(sys.argv[1])
