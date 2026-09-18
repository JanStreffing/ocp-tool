#!/bin/bash
#
# Automatic OASIS3-MCT installation script for ocp-tool
# Clones, builds, and installs OASIS3-MCT with pyOASIS support
#

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_DIR="${SCRIPT_DIR}/.oasis3-mct"
OASIS_REPO="https://git.smhi.se/jan.streffing/oasis3-mct-5"

echo "=========================================="
echo "OASIS3-MCT Installation for ocp-tool"
echo "=========================================="
echo ""
echo "Installation directory: ${INSTALL_DIR}"
echo "Repository: ${OASIS_REPO}"
echo ""

# Check prerequisites
echo "Checking prerequisites..."

# Check for conda environment
if [ -z "${CONDA_DEFAULT_ENV}" ]; then
    echo "Error: No conda environment active."
    echo "Please activate your conda environment first:"
    echo "  source ~/loadconda.sh"
    echo "  conda activate ocp-tool2"
    exit 1
fi

echo "✓ Conda environment: ${CONDA_DEFAULT_ENV}"

# Check for required compilers
if ! command -v mpifort &> /dev/null && ! command -v mpif90 &> /dev/null; then
    echo "Error: MPI Fortran compiler not found"
    echo "Install with: conda install -c conda-forge fortran-compiler openmpi"
    exit 1
fi

if ! command -v mpicc &> /dev/null; then
    echo "Error: MPI C compiler not found"
    echo "Install with: conda install -c conda-forge c-compiler openmpi"
    exit 1
fi

echo "✓ MPI compilers found"

# Check for NetCDF
if ! command -v nc-config &> /dev/null && ! command -v nf-config &> /dev/null; then
    echo "Error: NetCDF not found"
    echo "Install with: conda install -c conda-forge netcdf-fortran"
    exit 1
fi

echo "✓ NetCDF found"

# Clone OASIS3-MCT if not already present
if [ -d "${INSTALL_DIR}" ]; then
    echo ""
    read -p "OASIS3-MCT directory exists. Re-install? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing installation..."
        rm -rf "${INSTALL_DIR}"
    else
        echo "Using existing installation. Skipping to pyOASIS install..."
        SKIP_BUILD=1
    fi
fi

if [ -z "${SKIP_BUILD}" ]; then
    echo ""
    echo "Cloning OASIS3-MCT..."
    git clone "${OASIS_REPO}" "${INSTALL_DIR}"
    
    cd "${INSTALL_DIR}"
    
    echo ""
    echo "Configuring OASIS3-MCT with CMake..."
    
    # Detect compilers from conda environment
    MPIFORT=$(command -v mpifort || command -v mpif90)
    MPICC=$(command -v mpicc)
    MPICXX=$(command -v mpicxx || command -v mpic++)
    
    # Get NetCDF paths
    NETCDF_ROOT="${CONDA_PREFIX}"
    
    # Get Python paths for pyOASIS
    PYTHON_VERSION=$(python -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
    
    echo "Configuration:"
    echo "  MPI Fortran: ${MPIFORT}"
    echo "  MPI C:       ${MPICC}"
    echo "  NetCDF root: ${NETCDF_ROOT}"
    echo "  Python ver:  ${PYTHON_VERSION}"
    
    echo ""
    echo "Building OASIS3-MCT with CMake..."
    echo "This may take several minutes..."
    
    # Create build directory
    mkdir -p build
    cd build
    
    # Configure with CMake
    cmake .. \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
        -DCMAKE_Fortran_COMPILER="${MPIFORT}" \
        -DCMAKE_C_COMPILER="${MPICC}" \
        -DCMAKE_CXX_COMPILER="${MPICXX}" \
        -DNetCDF_ROOT="${NETCDF_ROOT}" \
        -DOASIS_USE_PYTHON=ON \
        -DCMAKE_BUILD_TYPE=Release
    
    if [ $? -ne 0 ]; then
        echo ""
        echo "Error: CMake configuration failed"
        exit 1
    fi
    
    # Build
    make -j4
    
    if [ $? -ne 0 ]; then
        echo ""
        echo "Error: OASIS3-MCT build failed"
        exit 1
    fi
    
    echo ""
    echo "✓ OASIS3-MCT built successfully"
fi

# Install pyOASIS
echo ""
echo "Installing pyOASIS..."

# Set required environment variables
export COUPLE="${INSTALL_DIR}"
export ARCH=conda
export OASIS_BUILD_DIR="${INSTALL_DIR}/build"

# pyOASIS is built as part of CMake if -DOASIS_USE_PYTHON=ON
# The shared library should be in the build directory
cd "${INSTALL_DIR}"

# Find the pyOASIS shared library
PYOASIS_LIB=$(find build -name "_pyoasis*.so" 2>/dev/null | head -1)

if [ -n "${PYOASIS_LIB}" ]; then
    echo "Found pyOASIS library: ${PYOASIS_LIB}"
    
    # Copy to Python site-packages
    SITE_PACKAGES=$(python -c "import sysconfig; print(sysconfig.get_path('purelib'))")
    cp "${PYOASIS_LIB}" "${SITE_PACKAGES}/"
    
    # Also copy the Python wrapper if it exists
    if [ -f "pyoasis/src/pyoasis.py" ]; then
        cp pyoasis/src/pyoasis.py "${SITE_PACKAGES}/"
    fi
    
    echo "✓ pyOASIS installed to ${SITE_PACKAGES}"
else
    echo "Warning: pyOASIS shared library not found"
    echo "Trying manual build..."
    
    cd pyoasis
    if [ -f "Makefile" ]; then
        make clean
        make
        
        # Find the built library
        PYOASIS_LIB=$(find . -name "_pyoasis*.so" 2>/dev/null | head -1)
        
        if [ -n "${PYOASIS_LIB}" ]; then
            SITE_PACKAGES=$(python -c "import sysconfig; print(sysconfig.get_path('purelib'))")
            cp "${PYOASIS_LIB}" "${SITE_PACKAGES}/"
            
            if [ -f "src/pyoasis.py" ]; then
                cp src/pyoasis.py "${SITE_PACKAGES}/"
            fi
            
            echo "✓ pyOASIS built and installed"
        else
            echo "Error: Failed to build pyOASIS"
            exit 1
        fi
    else
        echo "Error: pyoasis Makefile not found"
        exit 1
    fi
fi

# Verify installation
echo ""
echo "Verifying installation..."

python -c "import pyoasis; print('✓ pyOASIS import successful')" 2>/dev/null

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Installation Complete!"
    echo "=========================================="
    echo ""
    echo "OASIS3-MCT installed to: ${INSTALL_DIR}"
    echo "pyOASIS is now available in your conda environment"
    echo ""
    echo "You can now run ocp-tool with automatic RMP generation:"
    echo "  python run_ocp_tool.py configs/TCO95_CORE2.yaml"
    echo ""
else
    echo ""
    echo "Warning: pyOASIS import test failed"
    echo "You may need to set environment variables:"
    echo "  export COUPLE=${INSTALL_DIR}"
    echo "  export ARCH=conda"
fi
