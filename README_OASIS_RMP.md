# OASIS Remapping Weight Generation - Quick Guide

## What's New

OCP-Tool now **automatically generates OASIS remapping weight files** (`.rmp`) for OpenIFS-FESOM2 coupling!

### Features

✅ **Automatic FESOM grid export** - FESOM mesh written to OASIS format  
✅ **One-command OASIS installation** - Build OASIS3-MCT + pyOASIS automatically  
✅ **RMP file generation** - 14 coupling links (IFS↔FESOM) computed offline  
✅ **Optional feature** - Can be disabled if not needed  

## Quick Setup

### 1. Install OASIS3-MCT (One Time)

```bash
cd /work/ab0246/a270092/software/ocp-tool

# Activate conda environment
source ~/loadconda.sh
conda activate ocp-tool2

# Install build dependencies
conda install -c conda-forge fortran-compiler c-compiler cxx-compiler openmpi netcdf-fortran

# Run automatic installation
./install_oasis.sh
```

**Time:** ~5-10 minutes for full build

### 2. Run OCP-Tool

```bash
python run_ocp_tool.py configs/TCO95_CORE2.yaml
```

### Output Files

```
output/TCO95_CORE2/oasis_mct3_input/
├── grids.nc        # All grids: A095, L095, R095, feom, etc.
├── masks.nc        # Land-sea masks for all grids
├── areas.nc        # Cell areas for all grids
└── rmp_*.nc        # Remapping weight files (14 files)
```

## Without OASIS Installation

If you prefer not to install OASIS, the tool will still work:

1. **FESOM grid is still exported** to OASIS files
2. **RMP generation is skipped** with a warning
3. **Alternative workflows available**:
   - Let OASIS compute weights at runtime (slower)
   - Use ESMF_RegridWeightGen offline
   - Generate manually later

Disable in config:
```yaml
options:
  generate_rmp: false
```

## What Happens Under the Hood

### Step 4: Write OASIS Files
- OpenIFS grids: `A095`, `L095`, `R095`, `TCO95-land`
- Runoff grids: `RnfA`, `RnfO`
- **FESOM grid: `feom`** ← Uses pyfesom2 branch with OASIS export

### Step 4b: Generate RMP Files (Optional)
- Uses pyOASIS to compute remapping weights
- Creates `rmp_*.nc` for each coupling field
- Default: 10 IFS→FESOM + 4 FESOM→IFS fields

## Coupling Fields

**Atmosphere → Ocean** (10 fields):
- Radiative fluxes: SW, LW
- Heat fluxes: sensible, latent
- Freshwater: precipitation, snow, evaporation, sublimation
- Momentum: zonal/meridional wind stress

**Ocean → Atmosphere** (4 fields):
- SST (sea surface temperature)
- Ice fraction, temperature, albedo

## Requirements

### Always Required
- pyfesom2 (custom branch) - **Already installed** ✓
- NetCDF4, python-eccodes, etc. - From environment.yaml ✓

### For RMP Generation (Optional)
- OASIS3-MCT with pyOASIS - Install via `./install_oasis.sh`
- Fortran/C/C++ compilers - From conda-forge
- MPI - From conda-forge

## Troubleshooting

### "ModuleNotFoundError: No module named 'pyoasis'"

Run the installation script:
```bash
./install_oasis.sh
```

Or disable RMP generation in your config YAML.

### "FESOM grid not added"

Make sure you have the pyfesom2 branch installed:
```bash
cd /work/ab0246/a270092/software/pyfesom2
pip install -e .
```

### Build fails

Check you have compilers:
```bash
conda install -c conda-forge fortran-compiler c-compiler cxx-compiler openmpi netcdf-fortran
```

## Documentation

- **[Installation Guide](docs/INSTALL_OASIS.md)** - Detailed OASIS setup
- **[RMP Generation Guide](docs/OASIS_RMP_GENERATION.md)** - Full technical documentation
- **Installation script:** `install_oasis.sh`

## Implementation

**Modified files:**
- `ocp_tool/oasis_writer.py` - Added `_append_fesom_grid_to_oasis_files()`
- `ocp_tool/oasis_remap.py` - New module for pyOASIS wrapper
- `run_ocp_tool.py` - Added Step 4b (RMP generation)
- `ocp_tool/config.py` - Added `generate_rmp` option

**Uses pyfesom2 PR #223:**
- Branch: `/work/ab0246/a270092/software/pyfesom2`
- Function: `pyfesom2.write_fesom_oasis_files()`
- Converts FESOM mesh → OASIS grid format

**Uses custom OASIS3-MCT:**
- Repository: https://git.smhi.se/jan.streffing/oasis3-mct-5
- Installed to: `.oasis3-mct/` (inside ocp-tool)
- Provides: pyOASIS Python bindings

## Testing

Tested configuration: **TCO95_CORE2**
- ✓ FESOM grid export works
- ✓ OASIS files contain all grids including `feom`
- ✓ File sizes correct (+9.69 MB grids, +2.42 MB masks, +2.91 MB areas)
- ⚠ RMP generation requires pyOASIS installation

## Support

For issues or questions, see:
- `docs/INSTALL_OASIS.md` - Installation help
- `docs/OASIS_RMP_GENERATION.md` - Technical details
- GitHub PR: https://github.com/FESOM/pyfesom2/pull/223
