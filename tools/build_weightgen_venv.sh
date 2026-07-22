#!/bin/bash -l
# Build a Python venv whose mpi4py is compiled against the *module* OpenMPI 4.1.2,
# i.e. the same MPI that liboasis.cbind.so (pyoasis) is linked to. This is the
# interpreter the ocp-tool OASIS weight-generation worker (oasis_weights.run_worker)
# must run under. numpy/netCDF4 are inherited from the module python via
# --system-site-packages; only mpi4py is (re)built here.
set -euo pipefail

VENV=${1:-/work/ab0246/a270092/software/ocp-tool/env/weightgen-venv}
OASIS_BUILD_PATH=/work/ab0246/a270092/model_codes/awiesm3-develop-is/oasis

# Exact module set from comp-oasis3mct-5.0-smhi_script.sh
module purge
module load python3/2023.01-gcc-11.2.0
module load intel-oneapi-compilers/2022.0.1-gcc-11.2.0
module load intel-oneapi-mkl/2022.0.1-gcc-11.2.0
module load openmpi/4.1.2-intel-2021.5.0
module load hdf5/1.12.1-openmpi-4.1.2-intel-2021.5.0
module load netcdf-c/4.8.1-openmpi-4.1.2-intel-2021.5.0
module load netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0

echo "module python : $(which python3)"
echo "mpicc         : $(which mpicc)  ->  $(mpicc -show 2>/dev/null | head -c 60)"

# Fresh venv layered on the module python; keep numpy/netCDF4 from the module.
rm -rf "$VENV"
python3 -m venv --system-site-packages "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"

# Force a source build of mpi4py against the loaded OpenMPI mpicc. --ignore-installed
# so it lands in the venv even though the module python already ships an (MPICH) one;
# the venv site-packages precedes system site-packages on sys.path, so this wins.
# Override CC/LDSHARED to the OpenMPI wrappers so linking uses mpicc (which pulls
# libopen-pal/libopen-rte transitively) instead of the mambaforge compiler_compat/ld,
# which fails to resolve OpenMPI's internal opal_*/orte_* symbols.
export MPICC="$(which mpicc)"
export MPICXX="$(which mpicxx)"
export CC="$MPICC"
export CXX="$MPICXX"
export LDSHARED="$MPICC -shared"
echo "mpicc -show : $(mpicc -show 2>/dev/null | head -c 120)"
pip install --upgrade pip >/dev/null
# Latest mpi4py (new build backend, compatible with modern setuptools). Still a
# source build against $MPICC (OpenMPI 4.1.2); worker only uses `from mpi4py import MPI`.
pip install --no-binary=mpi4py --ignore-installed --no-deps mpi4py

echo "=== VERIFY: mpi4py MPI + co-load with pyoasis ==="
export OASIS_BUILD_PATH
export PYTHONPATH="$OASIS_BUILD_PATH/python:${PYTHONPATH:-}"
export MPIROOT="$(mpif90 -show 2>/dev/null | perl -lne 'm{ -I(.*?)/include } and print $1')"
export LD_LIBRARY_PATH="$MPIROOT/lib:$OASIS_BUILD_PATH/lib:${LD_LIBRARY_PATH:-}"
python3 - <<'PY'
from mpi4py import MPI
lib = MPI.Get_library_version().split(',')[0].strip()
print("mpi4py libver:", lib)
import pyoasis
print("pyoasis import OK:", pyoasis.__file__)
assert "Open MPI" in lib, "mpi4py is NOT OpenMPI -- build failed"
print("OK: mpi4py is OpenMPI and pyoasis co-loads")
PY
echo "venv ready: $VENV"
