#!/bin/bash -l
# Generate the natural-order OASIS remapping weights (rmp_*.nc) for TCO639 <-> NEXT.
# Drives the model's pyOASIS (one srun rank per coupling link) via ocp-tool's
# oasis_weights.generate_weights. Output rmp_*.nc land in OASIS_DIR.
#
#   feom in grids.nc / the produced rmp are in NATURAL (nod2d.out) node order.
#   Reorder to a partition (e.g. dist_9600) afterwards with submit_reorder_*.sh.
#
#SBATCH --job-name=ocp_weightgen_TCO639_NEXT
#SBATCH --account=ab0246
#SBATCH --partition=compute
#SBATCH --nodes=4                 # one whole node per link (4 links) -> concurrent
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

REPO=/work/ab0246/a270092/software/ocp-tool
VENV=$REPO/env/weightgen-venv
OASIS_BUILD_PATH=/work/ab0246/a270092/model_codes/awiesm3-develop-is/oasis
OASIS_DIR=$REPO/output/TCO639_NEXT/oasis_mct3_input
ATM_GRID=A640
THREADS=128

# --- Exact module env of the model's OASIS build (comp-oasis3mct-5.0-smhi) ------
module purge
module load python3/2023.01-gcc-11.2.0
module load intel-oneapi-compilers/2022.0.1-gcc-11.2.0
module load intel-oneapi-mkl/2022.0.1-gcc-11.2.0
module load openmpi/4.1.2-intel-2021.5.0
module load hdf5/1.12.1-openmpi-4.1.2-intel-2021.5.0
module load netcdf-c/4.8.1-openmpi-4.1.2-intel-2021.5.0
module load netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0

# venv with mpi4py built against this OpenMPI 4.1.2 (matches liboasis.cbind.so)
source "$VENV/bin/activate"

export OASIS_BUILD_PATH
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
export MPIROOT="$(mpif90 -show 2>/dev/null | perl -lne 'm{ -I(.*?)/include } and print $1')"
export LD_LIBRARY_PATH="$MPIROOT/lib:$OASIS_BUILD_PATH/lib:${LD_LIBRARY_PATH:-}"
# Big automatic arrays in OASIS/MCT (see oasis_weights notes); keep per-thread
# stack bounded so 128 threads don't consume excessive memory.
export OMP_STACKSIZE=64M

echo "=== ocp weightgen ${SLURM_JOB_ID} ==="
echo "host        : $(hostname)"
echo "nodes/links : ${SLURM_JOB_NUM_NODES} nodes"
echo "python      : $(which python)"
echo "OASIS_BUILD : ${OASIS_BUILD_PATH}"
echo "OASIS_DIR   : ${OASIS_DIR}"
echo "atm_grid    : ${ATM_GRID}   threads/link: ${THREADS}"
python -c "from mpi4py import MPI; print('mpi4py:', MPI.Get_library_version().split(',')[0].strip())"
echo

cd "$REPO"
python -c "
from ocp_tool.oasis_weights import generate_weights
files = generate_weights(
    '${OASIS_DIR}',
    atm_grid='${ATM_GRID}',
    method='existing',
    threads=${THREADS},
    oasis_build_path='${OASIS_BUILD_PATH}',
)
print('produced:', [f.name for f in files])
"

echo
echo "=== rmp files in OASIS_DIR ==="
ls -la "${OASIS_DIR}"/rmp_*.nc
