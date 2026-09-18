# Installing OASIS3-MCT for RMP Generation

## Quick Start

Automatic installation from Jan Streffing's OASIS3-MCT repository:

```bash
cd /work/ab0246/a270092/software/ocp-tool

# Activate conda environment
source ~/loadconda.sh
conda activate ocp-tool2

# Install build dependencies (if not already installed)
conda install -c conda-forge fortran-compiler c-compiler cxx-compiler openmpi netcdf-fortran

# Run automatic installation
./install_oasis.sh
```

The script will:
1. Clone OASIS3-MCT from https://git.smhi.se/jan.streffing/oasis3-mct-5
2. Configure for your conda environment
3. Build OASIS3-MCT libraries
4. Install pyOASIS Python bindings

**Installation directory:** `.oasis3-mct/` (inside ocp-tool)

## Build Dependencies

Required packages (installed via conda):

- **fortran-compiler** - Fortran compiler (gfortran)
- **c-compiler** - C compiler (gcc)
- **cxx-compiler** - C++ compiler (g++)
- **openmpi** - MPI implementation
- **netcdf-fortran** - NetCDF Fortran library

Install all at once:
```bash
conda install -c conda-forge fortran-compiler c-compiler cxx-compiler openmpi netcdf-fortran
```

## Verification

After installation, verify pyOASIS is available:

```bash
python -c "import pyoasis; print('pyOASIS available')"
```

## Using RMP Generation

Once OASIS is installed, ocp-tool will automatically generate remapping weight files:

```bash
python run_ocp_tool.py configs/TCO95_CORE2.yaml
```

Output:
```
Step 4b: Generating OASIS remapping weight files...
  Initializing OASIS component...
  Processing coupling link 1/14: A095 -> feom (SISUTESW)
  ...
  ✓ Generated 14 remapping weight files
```

## Troubleshooting

### Build Fails with Compiler Errors

Make sure you've activated the conda environment:
```bash
source ~/loadconda.sh
conda activate ocp-tool2
```

Check compiler availability:
```bash
which mpifort  # or mpif90
which mpicc
```

### pyOASIS Import Fails

Set environment variables manually:
```bash
export COUPLE=/work/ab0246/a270092/software/ocp-tool/.oasis3-mct
export ARCH=conda
```

### Git Clone Fails

If you don't have access to git.smhi.se, contact Jan Streffing for repository access or use the official OASIS3-MCT from CERFACS:
```bash
# Edit install_oasis.sh and change:
OASIS_REPO="https://oasis.cerfacs.fr/git/oasis3-mct.git"
```

## Manual Installation

If the automatic script fails, you can build manually:

```bash
# Clone repository
git clone https://git.smhi.se/jan.streffing/oasis3-mct-5 .oasis3-mct
cd .oasis3-mct

# Set environment
export COUPLE=$(pwd)
export ARCH=conda

# Edit util/make_dir/make.conda (created by script)
# Then build:
make -f TopMakefileOasis3 realclean
make -f TopMakefileOasis3

# Install pyOASIS
cd pyoasis
python setup.py install
```

## Disabling RMP Generation

If you prefer not to install OASIS, set in your config YAML:

```yaml
options:
  generate_rmp: false  # Skip automatic RMP generation
```

The FESOM grid will still be exported to OASIS files, and you can generate RMP files later using alternative methods (ESMF_RegridWeightGen or runtime OASIS).
