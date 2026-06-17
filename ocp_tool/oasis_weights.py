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
    """One OASIS coupling link == one remapping weight file rmp_<src>_to_<tgt>_<MAP>.nc.

    For weight generation only the single ``SCRIPR`` remapping transform matters
    (LOCTRANS time transforms and post-remap CONSERV conservation do not create
    rmp files), so a link carries just its SCRIPR method line, e.g.
    ``"GAUSWGT D SCALAR LATITUDE 1 25 0.1"`` or ``"BICUBIC D SCALAR LATITUDE 15"``.
    """

    source: str   # source grid name in grids.nc (e.g. "A096")
    target: str   # target grid name (e.g. "feom")
    scrip: str     # SCRIPR method line (map kind + its parameters)

    @property
    def map_name(self) -> str:
        return self.scrip.split()[0]  # GAUSWGT / BICUBIC / CONSERV / DISTWGT

    def as_dict(self):
        return {"source": self.source, "target": self.target, "scrip": self.scrip}

    @staticmethod
    def from_dict(d):
        # tolerate the older {"scripr": [...]} form
        if "scrip" in d:
            return Link(d["source"], d["target"], d["scrip"])
        scripr = d["scripr"]
        method = next((s for s in scripr if s.split()[0] != "SCRIPR"), scripr[-1])
        return Link(d["source"], d["target"], method)


# Default awiesm3 feom links -- these reproduce the remapping methods of the
# existing TCO95/CORE2 namcouple (non-conservative, centre based). When the
# feom grid is regenerated with corners, set ``method="conserv"`` to switch the
# flux links to conservative remapping.
def awiesm3_feom_links(atm_grid: str = "A096", method: str = "existing") -> List[Link]:
    if method not in ("existing", "conserv"):
        raise ValueError(f"unknown method profile: {method}")
    gauswgt = "GAUSWGT D SCALAR LATITUDE 1 25 0.1"
    gauswgt_u = "GAUSWGT U SCALAR LATITUDE 1 25 0.1"
    bicubic = "BICUBIC D SCALAR LATITUDE 15"
    gauswgt_lr = "GAUSWGT LR SCALAR LATITUDE 1 25 0.1"
    conserv_lr = "CONSERV LR SCALAR LATITUDE 1 25 0.1"

    flux_atm_to_oce = conserv_lr if method == "conserv" else gauswgt
    rnf_to_oce = conserv_lr if method == "conserv" else gauswgt_lr

    return [
        Link(atm_grid, "feom", flux_atm_to_oce),   # A096 -> feom (heat/freshwater fluxes)
        Link(atm_grid, "feom", bicubic),           # A096 -> feom (state, bicubic)
        Link("feom", atm_grid, gauswgt_u),         # feom -> A096 (SST, ice state)
        Link("RnfO", "feom", rnf_to_oce),          # runoff -> feom
    ]


# ---------------------------------------------------------------------------
# Grid geometry read back from the OASIS description files
# ---------------------------------------------------------------------------
def _copy_without_grid(src, dst, prefix: str = "feom"):
    """Copy an OASIS grid-description file, dropping one grid's vars and dims.

    Used to build a feom-free template before pyfesom2 writes the new feom grid:
    pyfesom2's incremental writer refuses to change an existing ``x_feom`` size,
    so the old feom (vars ``feom.*`` and dims ``*_feom``) must be removed first.
    Atmosphere (A096) and runoff (RnfO) grids are preserved.
    """
    with Dataset(src) as ds_in, Dataset(dst, "w", format=ds_in.file_format) as ds_out:
        ds_out.setncatts({k: ds_in.getncattr(k) for k in ds_in.ncattrs()})
        drop_dims = {d for d in ds_in.dimensions if d.endswith(f"_{prefix}")}
        for name, dim in ds_in.dimensions.items():
            if name in drop_dims:
                continue
            ds_out.createDimension(name, None if dim.isunlimited() else len(dim))
        for name, var in ds_in.variables.items():
            if name.startswith(f"{prefix}.") or drop_dims.intersection(var.dimensions):
                continue
            out = ds_out.createVariable(name, var.datatype, var.dimensions)
            out.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
            out[:] = var[:]


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
def component_names(nlinks: int, prefix: str = "ocpw") -> List[str]:
    """OASIS component names, one per link/rank (must match run_worker)."""
    return [f"{prefix}{i:02d}" for i in range(nlinks)]


def write_namcouple(run_dir: Path, links: List[Link], runtime: int = 3600):
    """Write a namcouple that makes OASIS compute one SCRIPR weight file per link.

    Format mirrors the model's own namcouple exactly: $NBMODEL with the per-rank
    component names, then for each field the grid-size line
    ``<snx> <sny> <tnx> <tny> <src> <tgt> LAG=0`` and a single SCRIPR transform.
    """
    comps = component_names(len(links))
    lines = [
        "# namcouple generated by ocp_tool.oasis_weights for weight generation only",
        " $NFIELDS", f"   {len(links)}", " $END",
        " $NBMODEL", f"   {len(comps)} " + " ".join(comps), " $END",
        " $RUNTIME", f"   {runtime}", " $END",
        " $NLOGPRT", "   1 0 0", " $END",
        " $NNOREST", "   T", " $END",
        " $STRINGS",
    ]
    for i, lk in enumerate(links):
        src = OasisFileGrid(lk.source, run_dir)
        tgt = OasisFileGrid(lk.target, run_dir)
        snx, sny = src.shape
        tnx, tny = tgt.shape
        lines += [
            "#",
            f"# link {i}: {lk.source} -> {lk.target} ({lk.map_name})",
            f"VAR_{i:02d}_S VAR_{i:02d}_T 1 {runtime} 1 none EXPORTED",
            f"{snx} {sny} {tnx} {tny} {lk.source} {lk.target} LAG=0",
            "P 0 P 0",
            "SCRIPR",
            lk.scrip,
        ]
    lines += ["#", " $END"]
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
    threads: int = 8,
    launcher: Optional[List[str]] = None,
    worker_python: Optional[str] = None,
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

    # OASIS computes the SCRIP weights with OpenMP; the per-link cost is
    # dominated by the GAUSWGT nearest-neighbour search over the ~10^5 ocean
    # nodes, which threads well. Give each rank `threads` cores.
    threads = max(1, int(threads))
    env["OASIS_OMP_NUM_THREADS"] = str(threads)
    env["OMP_NUM_THREADS"] = str(threads)

    if launcher is None:
        if subprocess.call(["bash", "-lc", "command -v srun >/dev/null 2>&1"]) == 0:
            launcher = ["srun", "-n", str(nlinks),
                        f"--cpus-per-task={threads}", "--cpu-bind=cores"]
        else:
            launcher = ["mpirun", "-n", str(nlinks)]

    # The srun worker only needs pyoasis + (OpenMPI) mpi4py + netCDF4 -- it reads
    # grids from files, no pyfesom2 -- so it may run a different interpreter than
    # this driver (which needs pyfesom2 to regenerate the feom grid).
    py = worker_python or sys.executable
    cmd = launcher + [py, "-m", "ocp_tool.oasis_weights", str(oasis_dir)]
    print(f"[oasis_weights] launching: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=oasis_dir, env=env)

    produced = sorted(p.name for p in oasis_dir.glob("rmp_*.nc"))
    print(f"[oasis_weights] produced {len(produced)} weight files: {produced}")
    return [oasis_dir / n for n in produced]


# ---------------------------------------------------------------------------
# High-level entry point for dynamic-mesh coupling (called from ice2fesom)
# ---------------------------------------------------------------------------
def regenerate_for_mesh(
    mesh_path,
    template_oasis_dir,
    out_dir,
    *,
    method: str = "existing",
    atm_grid: str = "A096",
    threads: int = 8,
    worker_python: Optional[str] = None,
    oasis_build_path: Optional[str] = None,
):
    """Rebuild the OASIS grid files + feom remap weights for a new FESOM mesh.

    The driver (this call) needs pyfesom2 to write the feom grid; the srun
    weight-gen worker needs (OpenMPI) mpi4py + pyoasis. If those live in
    different interpreters, pass ``worker_python`` for the latter.

    Steps:
      1. copy grids.nc/masks.nc/areas.nc from ``template_oasis_dir`` (which holds
         the atmosphere A096 and runoff RnfO grids) into ``out_dir``;
      2. overwrite the ``feom`` grid in them for ``mesh_path`` (the new submesh),
         *with corners*, via pyfesom2.write_fesom_oasis_files;
      3. run pyOASIS to (re)compute rmp_*feom*.nc for the new node count.

    Returns the list of produced rmp_*.nc paths. ``out_dir`` then contains a
    consistent grids/masks/areas + rmp_* set to stage for the next FESOM run.
    """
    template = Path(template_oasis_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Copy the atmosphere/runoff grids but drop any stale feom (its dims would
    # otherwise block pyfesom2 from writing the new feom node count).
    for f in ("grids.nc", "masks.nc", "areas.nc"):
        _copy_without_grid(template / f, out_dir / f, prefix="feom")

    # Overwrite the feom grid (centres + 4 corners + area + mask) for the new mesh.
    import pyfesom2 as pf

    mesh = pf.load_mesh(str(Path(mesh_path)))
    pf.write_fesom_oasis_files(mesh=mesh, output_dir=str(out_dir),
                               prefix="feom", overwrite=True)

    return generate_weights(
        out_dir, atm_grid=atm_grid, method=method, threads=threads,
        worker_python=worker_python, oasis_build_path=oasis_build_path,
    )


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("usage: python -m ocp_tool.oasis_weights <run_dir>\n")
        sys.exit(2)
    run_worker(sys.argv[1])
