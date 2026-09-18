# OASIS Remapping Weight Generation for OCP-Tool

## Overview

OCP-Tool now supports automatic generation of OASIS remapping weight files (`.rmp` files) for coupling OpenIFS and FESOM2. This feature uses pyOASIS (OASIS3-MCT Python bindings) to compute conservative remapping weights offline.

## What's New

### Automatic FESOM Grid Export

When processing configurations with FESOM ocean grids (e.g., CORE2, CORE3, DART), ocp-tool now automatically:

1. **Exports FESOM mesh to OASIS format** after writing OpenIFS grids
2. **Adds FESOM grid** (`feom`) to `grids.nc`, `masks.nc`, `areas.nc`
3. **Uses pyfesom2** branch with `write_fesom_oasis_files()` functionality

**Grids in OASIS files:**
- `A###` - Atmosphere grid (cell centers)
- `L###` - Land grid (same as atmosphere)
- `R###` - Runoff grid
- `TCO##-land` - LPJ-GUESS vegetation grid
- `RnfA`, `RnfO` - Runoff mapper grids
- `feom` - FESOM ocean grid ✓ **NEW**

### Optional RMP File Generation

If pyOASIS is available, ocp-tool can automatically generate remapping weight files for all coupling links between OpenIFS and FESOM.

## Requirements

### pyfesom2 with OASIS Export Support

**Already installed** via the custom branch at `/work/ab0246/a270092/software/pyfesom2`:

```bash
cd /work/ab0246/a270092/software/pyfesom2
conda activate ocp-tool2
pip install -e .
```

This provides the `write_fesom_oasis_files()` function used by ocp-tool.

### pyOASIS (Optional - For RMP Generation)

pyOASIS requires OASIS3-MCT to be compiled with Python bindings enabled.

#### Installation

1. **Download OASIS3-MCT** (version ≥ 5.0)
   ```bash
   wget https://oasis.cerfacs.fr/download/oasis3-mct_5.0.tar.gz
   tar xzf oasis3-mct_5.0.tar.gz
   cd oasis3-mct_5.0
   ```

2. **Configure with Python bindings**
   ```bash
   # Edit util/make_dir/make.inc
   # Set: USE_PYTHON = yes
   # Set: PYTHON_INC = -I/path/to/python/include
   # Set: PYTHON_LIB = -L/path/to/python/lib -lpython3.x
   ```

3. **Compile OASIS3-MCT**
   ```bash
   make realclean
   make
   ```

4. **Install pyOASIS**
   ```bash
   cd pyoasis
   python setup.py install
   ```

#### Verification
```bash
python -c "import pyoasis; print('pyOASIS available')"
```

## Usage

### Configuration

Enable/disable RMP generation in your config YAML:

```yaml
options:
  verbose: true
  parallel_workers: 4
  use_dask: true
  generate_rmp: true  # Set to false to skip RMP generation
```

### Running OCP-Tool

```bash
cd /work/ab0246/a270092/software/ocp-tool
conda activate ocp-tool2
python run_ocp_tool.py configs/TCO95_CORE2.yaml
```

#### Processing Steps

**Step 4:** Write OASIS grid/mask/area files
- OpenIFS grids (A095, L095, R095, TCO95-land)
- Runoff mapper grids (RnfA, RnfO)
- **FESOM grid (feom)** ✓

**Step 4b:** Generate OASIS remapping weights
- Creates `rmp_*.nc` files for each coupling link
- Default: 14 coupling links (10 IFS→FESOM, 4 FESOM→IFS)
- Requires pyOASIS installation

### Output Files

After successful processing:

```
output/TCO95_CORE2/oasis_mct3_input/
├── grids.nc        # Grid coordinates + corners (all grids including feom)
├── masks.nc        # Land-sea masks (all grids including feom)
├── areas.nc        # Cell areas (all grids including feom)
├── rmp_*.nc        # Remapping weight files (if pyOASIS available)
└── namcouple       # OASIS coupling configuration (if generated)
```

## Default Coupling Links

### Atmosphere → Ocean (10 fields)

| Field | Description | Remapping |
|-------|-------------|-----------|
| SISUTESW | Shortwave radiation | DISTWGT |
| SISUTELW | Longwave radiation | DISTWGT |
| SISUTESENS | Sensible heat flux | DISTWGT |
| SISUTELAT | Latent heat flux | DISTWGT |
| SISUTEPRE | Precipitation | DISTWGT |
| SISUTESNOW | Snow | DISTWGT |
| SISUTESUB | Sublimation | DISTWGT |
| SISUTEEVAP | Evaporation | DISTWGT |
| SISUTAUXI | Zonal wind stress | DISTWGT |
| SISUTAUYI | Meridional wind stress | DISTWGT |

### Ocean → Atmosphere (4 fields)

| Field | Description | Remapping |
|-------|-------------|-----------|
| A_SST | Sea surface temperature | DISTWGT |
| A_Ice_frac | Sea ice fraction | DISTWGT |
| A_Ice_temp | Sea ice temperature | DISTWGT |
| A_Ice_albedo | Sea ice albedo | DISTWGT |

## Alternative Workflows (Without pyOASIS)

### Option 1: Runtime Weight Computation

Let OASIS3-MCT compute weights during the first coupling timestep:
- Set `generate_rmp: false` in config
- OASIS will create `.rmp` files automatically on first run
- **Slower** but no offline setup needed

### Option 2: ESMF_RegridWeightGen

Use ESMF's offline weight generation tool (faster than OASIS):

```bash
# Install ESMF
conda install -c conda-forge esmf

# Generate weights
ESMF_RegridWeightGen \
  --source grids.nc \
  --destination grids.nc \
  --src_var A095.lon,A095.lat \
  --dst_var feom.lon,feom.lat \
  --method conserve \
  --weight rmp_A095_to_feom_conserve.nc
```

### Option 3: Manual pyOASIS Script

Create a standalone Python script using pyOASIS:

```python
import pyoasis
from pathlib import Path
import os

oasis_dir = Path("output/TCO95_CORE2/oasis_mct3_input")
os.chdir(oasis_dir)

# Initialize OASIS component
comp = pyoasis.Component("weight_gen")

# Define partitions (serial mode)
partition_ifs = pyoasis.SerialPartition(0)
partition_fesom = pyoasis.SerialPartition(0)

# Define coupling variable
var_src = pyoasis.Var("SST_src", partition_ifs, bundle_size=1, io_type=pyoasis.OasisOut)
var_tgt = pyoasis.Var("SST_tgt", partition_fesom, bundle_size=1, io_type=pyoasis.OasisIn)

# Trigger weight computation
comp.enddef()

# Terminate
comp.terminate()
```

## Troubleshooting

### FESOM Grid Not Added

**Error:** `Warning: pyfesom2 not available`

**Solution:**
```bash
cd /work/ab0246/a270092/software/pyfesom2
conda activate ocp-tool2
pip install -e .
```

### pyOASIS Import Error

**Error:** `ModuleNotFoundError: No module named 'pyoasis'`

**Solution:** Either:
1. Install pyOASIS (see Requirements section)
2. Set `generate_rmp: false` in config
3. Use alternative workflow (ESMF or runtime generation)

### MPI Errors During RMP Generation

pyOASIS may require MPI even for serial execution:

```bash
# Run with mpirun if needed
mpirun -n 1 python run_ocp_tool.py configs/TCO95_CORE2.yaml
```

## Implementation Details

### Code Structure

- `ocp_tool/oasis_writer.py` - OASIS grid file generation + FESOM grid export
- `ocp_tool/oasis_remap.py` - pyOASIS wrapper for RMP generation
- `run_ocp_tool.py` - Main pipeline with Step 4b (RMP generation)

### FESOM Grid Export

Uses pyfesom2's `write_fesom_oasis_files()`:
- Loads FESOM mesh via `pyfesom2.load_mesh()`
- Exports node coordinates as grid centers
- Exports element centroids as grid corners (4 per node)
- Computes node areas by distributing element areas
- All nodes marked as ocean (mask=1)

### Grid Prefix Convention

- OpenIFS: `A###`, `L###`, `R###`
- FESOM: `feom`
- Runoff: `RnfA`, `RnfO`

## Testing

Tested with:
- Configuration: `TCO95_CORE2.yaml`
- FESOM mesh: CORE2 (126,858 nodes, 244,659 elements)
- pyfesom2 branch: `/work/ab0246/a270092/software/pyfesom2`

**Results:**
- ✓ FESOM grid successfully added to OASIS files
- ✓ File sizes: grids.nc +9.69 MB, masks.nc +2.42 MB, areas.nc +2.91 MB
- ✓ Grid `feom` appears in all three OASIS files
- ⚠ RMP generation requires pyOASIS installation

## References

- [OASIS3-MCT Documentation](https://oasis.cerfacs.fr/en/)
- [pyfesom2 PR #223](https://github.com/FESOM/pyfesom2/pull/223) - OASIS export functionality
- [rdy2cpl Tool](https://github.com/FESOM/rdy2cpl) - Original inspiration for RMP generation
- [ESMF RegridWeightGen](https://www.earthsystemcog.org/projects/regridweightgen/)
