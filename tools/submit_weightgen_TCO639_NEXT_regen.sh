#!/bin/bash -l
# Regenerate ONLY the two rmp links that ocp-tool's old defaults got wrong for
# the awiesm3 (N56-style) namcouple:
#   feom -> A640  GAUSWGT U 1 9 2   (was 25 nn; SST/ice/currents back to atm)
#   R640 -> RnfA  GAUSWGT D 1 25 0.1 (was missing; runoff/calving to mapper)
# The other three (A640->feom GAUSWGT/BICUBIC, RnfO->feom GAUSWGT) are already
# correct and are left in place. Output rmp land in OASIS_DIR alongside them.
#
#SBATCH --job-name=ocp_weightgen_regen_TCO639_NEXT
#SBATCH --account=ab0246
#SBATCH --partition=compute
#SBATCH --nodes=2                 # one node per link (2 links) -> concurrent
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

REPO=/work/ab0246/a270092/software/ocp-tool
VENV=$REPO/env/weightgen-venv
OASIS_BUILD_PATH=/work/ab0246/a270092/model_codes/awiesm3-develop-is/oasis
OASIS_DIR=$REPO/output/TCO639_NEXT/oasis_mct3_input
THREADS=128

module purge
module load python3/2023.01-gcc-11.2.0
module load intel-oneapi-compilers/2022.0.1-gcc-11.2.0
module load intel-oneapi-mkl/2022.0.1-gcc-11.2.0
module load openmpi/4.1.2-intel-2021.5.0
module load hdf5/1.12.1-openmpi-4.1.2-intel-2021.5.0
module load netcdf-c/4.8.1-openmpi-4.1.2-intel-2021.5.0
module load netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0
source "$VENV/bin/activate"

export OASIS_BUILD_PATH
export PYTHONPATH="$REPO:${PYTHONPATH:-}"
export MPIROOT="$(mpif90 -show 2>/dev/null | perl -lne 'm{ -I(.*?)/include } and print $1')"
export LD_LIBRARY_PATH="$MPIROOT/lib:$OASIS_BUILD_PATH/lib:${LD_LIBRARY_PATH:-}"
export OMP_STACKSIZE=64M

echo "=== ocp weightgen regen ${SLURM_JOB_ID} on $(hostname) ==="
python -c "from mpi4py import MPI; print('mpi4py:', MPI.Get_library_version().split(',')[0].strip())"

cd "$REPO"
python -c "
from ocp_tool.oasis_weights import generate_weights, Link
links = [
    Link('feom', 'A640', 'GAUSWGT U SCALAR LATITUDE 1 9 2'),
    Link('R640', 'RnfA', 'GAUSWGT D SCALAR LATITUDE 1 25 0.1'),
]
files = generate_weights('${OASIS_DIR}', links=links, threads=${THREADS},
                         oasis_build_path='${OASIS_BUILD_PATH}')
print('produced:', [f.name for f in files])
"

echo
echo "=== full rmp set in OASIS_DIR now ==="
ls -la "${OASIS_DIR}"/rmp_*.nc
