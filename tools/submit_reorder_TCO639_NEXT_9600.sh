#!/bin/bash
# Reorder the natural-order rmp_*.nc for TCO639<->NEXT into dist_9600 (9600-rank)
# FESOM partition order, using the prebuilt reorder_oasis tool.
#
#   DIST_OLD = dist_natural  (identity rpart.out: mapping(n)=n, i.e. nod2d order)
#   DIST_NEW = dist_9600     (the partition the coupled run will use)
#
# Only rmp_* / rstos* / rstas* / vegin are touched (feom-side addresses permuted);
# grids/masks/areas are NOT reordered by this tool (feom there is added/ordered at
# runtime for the chosen partition).
#
#SBATCH --job-name=ocp_reorder_TCO639_NEXT_9600
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

REORDER_DIR=/work/ab0246/a270092/input/oasis/cy48r1/oasis_reorder_tool
MESH=/work/ab0246/a270092/input/fesom2/next/mesh
DIST_OLD=$MESH/dist_natural
DIST_NEW=$MESH/dist_9600
OASIS_IN=/work/ab0246/a270092/software/ocp-tool/output/TCO639_NEXT/oasis_mct3_input
OASIS_OUT=/work/ab0246/a270092/software/ocp-tool/output/TCO639_NEXT/oasis_mct3_input_dist9600

module purge
module load intel-oneapi-compilers/2023.2.1-gcc-11.2.0
module load intel-oneapi-mpi/2021.5.0-intel-2021.5.0
module load netcdf-fortran/4.5.3-intel-oneapi-mpi-2021.5.0-intel-2021.5.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OMP_PROC_BIND=close
export OMP_PLACES=cores
ulimit -s unlimited

mkdir -p "$OASIS_OUT"

echo "=== reorder TCO639_NEXT natural -> dist_9600 (${SLURM_JOB_ID}) ==="
echo "DIST_OLD : $DIST_OLD"
echo "DIST_NEW : $DIST_NEW"
echo "OASIS_IN : $OASIS_IN"
echo "OASIS_OUT: $OASIS_OUT"
echo "binary   : $REORDER_DIR/reorder_oasis"
echo

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
srun --mpi=pmi2 --chdir="$REORDER_DIR" ./reorder_oasis "$DIST_OLD" "$DIST_NEW" "$OASIS_IN" "$OASIS_OUT"

echo
echo "=== output ==="
ls -la "$OASIS_OUT/"
