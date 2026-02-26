# Paleo Modifications Plan for OCP-Tool

Replaces the 5 bash scripts (`mod_input_oifs.sh`, `mod_topo_oifs.sh`, `convert_topodata.sh`, `create_subgrid_oro.sh`, `mod_oro_oifs.sh`) and the calnoro Fortran workflow with pure Python inside ocp-tool. No CDO, no bash.

## Overview of Bash Script → Python Mapping

| Bash Script | Python Module | Purpose |
|---|---|---|
| `mod_input_oifs.sh` | `ocp_tool/paleo_input.py` | Modify ICMGG land surface variables from reconstructions |
| `mod_topo_oifs.sh` | `ocp_tool/paleo_topo.py` | Modify ICMSH topography (anomaly method + spectral) |
| `convert_topodata.sh` + calnoro + `create_subgrid_oro.sh` | `ocp_tool/paleo_subgrid_oro.py` | Compute subgrid-scale orography via calnoro |
| `mod_oro_oifs.sh` | `ocp_tool/paleo_subgrid_oro.py` | Apply SSO parameters to ICMGG |

## Libraries

- **eccodes / gribapi** — GRIB I/O (already in ocp-tool)
- **scipy.interpolate** — replace `cdo remapycon` (conservative) / `remapdis` (distance-weighted)
- **numpy** — array operations, value remapping (replace `cdo setvals`, `setrtoc`, masking)
- **xarray + cfgrib** — convenient GRIB ↔ NetCDF reading (optional, for reconstruction NetCDF files)
- **subprocess** — call compiled calnoro binary

---

## Checklist

### 1. Config & Lookup Tables
- [x] Extend `OCPConfig` with `PaleoConfig` dataclass: paths to reconstruction files (ice mask, topography, LSM, lake, vegetation/biome, soil), path to modern reference files, calnoro binary path, experiment ID (`lsm_id`)
- [x] Define lookup tables as Python dicts:
  - PlioMIP3 soil type → OpenIFS soil category (code 43)
  - OpenIFS soil category → volumetric soil water per layer (codes 39–42)
  - PlioMIP3 biome → OpenIFS high vegetation type (code 30)
  - PlioMIP3 biome → OpenIFS low vegetation type (code 29)
  - High vegetation type → cover fraction (code 28)
  - Low vegetation type → cover fraction (code 27)
  - High vegetation type → LAI (code 67)
  - Low vegetation type → LAI (code 66)

### 2. Mask Creation (`paleo_input.py` — Part 1)
Replaces `mod_input_oifs.sh` lines 58–186.

**No intermediate grid needed.** The polygon method (`read_fesom_grid_polygon`) already produces the paleo LSM directly on the OpenIFS Gaussian grid. Derivative masks are computed directly via numpy on the Gaussian grid, and reconstruction data is interpolated directly to the Gaussian grid via `scipy.interpolate.griddata`.

- [x] Read pre-modification (CORE2) and post-modification (paleo) LSM from GRIB on Gaussian grid
- [x] Create ocean/land masks directly: `land = paleo_lsm >= 0.5`, `ocean = paleo_lsm < 0.5`
- [x] Create derivative masks via numpy:
  - `added_land = (core2_lsm < 0.5) & (paleo_lsm >= 0.5)` — ocean→land
  - `drowned = (core2_lsm >= 0.5) & (paleo_lsm < 0.5)` — land→ocean
- [x] Read reconstruction ice mask, interpolate directly to Gaussian grid (`scipy.interpolate.griddata`), create "reduced ice sheet" mask by comparing with modern ice field

### 3. Ice / Snow Depth Modification (`paleo_input.py` — Part 2)
Replaces `mod_input_oifs.sh` lines 187–202.
- [x] Read reconstruction ice mask
- [x] Map ice values: `1→0` (no ice), `2→10` (ice sheet = snow depth 10 m water equiv.)
- [x] Remap to OpenIFS grid, assign GRIB code 141

### 4. Lake Modification (`paleo_input.py` — Part 3)
Replaces `mod_input_oifs.sh` lines 204–242.
- [x] Read modern + paleo lake masks (value 1100 = lake fraction in PlioMIP3)
- [x] Compute lake anomaly, add to existing OpenIFS lake field (code 26)
- [x] Clamp to binary 0/1, remap to OpenIFS grid

### 5. Soil Modification (`paleo_input.py` — Part 4)
Replaces `mod_input_oifs.sh` lines 244–301.
- [x] Read reconstruction soil data
- [x] Distance-weighted fill relative to land mask
- [x] Apply PlioMIP3 → OpenIFS soil category lookup
- [x] Remap to OpenIFS grid, assign GRIB code 43

### 6. Volumetric Soil Water Layers (`paleo_input.py` — Part 5)
Replaces `mod_input_oifs.sh` lines 303–342.
- [x] From soil categories, apply per-layer lookup tables:
  - Layer 1 (code 39, level 0): `{1: 0.210, 2: 0.283, 3: 0.373, 4: 0.209}`
  - Layer 2 (code 40, level 7): `{1: 0.212, 2: 0.279, 3: 0.374, 4: 0.248}`
  - Layer 3 (code 41, level 28): `{1: 0.208, 2: 0.270, 3: 0.375, 4: 0.262}`
  - Layer 4 (code 42, level 100): `{1: 0.188, 2: 0.300, 3: 0.399, 4: 0.255}`

### 7. Vegetation Type Modification (`paleo_input.py` — Part 6)
Replaces `mod_input_oifs.sh` lines 344–406.
- [x] Read reconstruction biome data
- [x] Distance-weighted fill relative to land mask
- [x] Apply biome → high vegetation type lookup: `{1→6, 2→3, 6→5, 7→4}`
- [x] Apply biome → low vegetation type lookup: `{3→7, 4→2, 5→11, 8→9}`
- [x] Remap both to OpenIFS grid, assign codes 30 and 29

### 8. Vegetation Cover Modification (`paleo_input.py` — Part 7)
Replaces `mod_input_oifs.sh` lines 408–444.
- [x] High veg cover from type: `{3→0.95, 4→0.925, 5→0.95, 6→0.95}` → code 28
- [x] Low veg cover from type: `{2→0.87, 7→0.90, 9→0.40, 11→0.10}` → code 27

### 9. Vegetation LAI Modification (`paleo_input.py` — Part 8)
Replaces `mod_input_oifs.sh` lines 446–482.
- [x] High veg LAI from type: `{3→5.9, 4→4.89, 5→5.97, 6→6.44}` → code 67
- [x] Low veg LAI from type: `{2→2.92, 7→3.79, 9→2.95, 11→3.11}` → code 66

### 10. Remaining Fields Interpolation (`paleo_input.py` — Part 9)
Replaces `mod_input_oifs.sh` lines 484–608.
- [x] For each remaining GRIB code (139, 170, 183, 236, 198, 235, 31, 8–14, 238, 32–38, 148, 174, 15–18, 74):
  - Mask to ocean-only, distance-weighted fill for drowned points
  - Mask to land-only, distance-weighted fill for added land + reduced ice sheet
  - Combine ocean + land fields
  - Remap to OpenIFS grid
- [x] Merge all modified fields into single GRIB
- [x] Multiply by template ("ones") to fix codes/levels
- [x] Replace fields in ICMGGaackINIT → produces `ICMGGaackINIT_<lsm_id>_MASKS`

### 11. Topography Modification (`paleo_topo.py`)
Replaces `mod_topo_oifs.sh`.
- [x] Read modern + paleo topography from reconstruction files
- [x] Clip below sea level to 0
- [x] Distance-weighted fill for land mask
- [x] Compute topography anomaly (paleo − modern)
- [x] Extract existing OpenIFS spectral topography from ICMSHaackINIT (sp2gp, ÷9.8)
- [x] Add anomaly, mask by land, clip below 0
- [x] Convert back to spectral (gp2sp), ×9.8, remap to TCO grid
- [x] Replace in ICMSHaackINIT → produces `ICMSHaackINIT_<lsm_id>`

### 12. Subgrid-Scale Orography via Calnoro (`paleo_subgrid_oro.py`)
Replaces `convert_topodata.sh` + calnoro + `create_subgrid_oro.sh` + `mod_oro_oifs.sh`.
- [x] **Prepare calnoro input**: remap paleo topography to r720x360, rename variable to OROMEA, write as `topodata.nc`
- [x] **Run calnoro**: call compiled binary via `subprocess`, feed resolution (63)
- [x] **Post-process output**: read `sso_par_fil.srv` (SERVICE format), set T63 grid, flip latitudes, convert to NetCDF, rename variables to OROMEA/OROSTD/OROSIG/OROGAM/OROTHE/OROPIC/OROVAL
- [x] **Apply to ICMGG**: for each SSO variable:
  - OROSTD → code 160 (standard deviation of orography)
  - OROSIG → code 163 (slope, ×10)
  - OROGAM → code 161 (anisotropy)
  - OROTHE → code 162 (angle, ×π/180)
- [x] Remap to OpenIFS grid, merge, replace in `ICMGGaackINIT_<lsm_id>_MASKS` → final `ICMGGaackINIT_<lsm_id>_MASKS_ORO`

### 13. Integration
- [x] Add paleo steps to `run_ocp_tool.py` (conditional on `paleo.enabled` in config)
- [x] Add example paleo config YAML (e.g., `configs/TCO95_CORE2_EP.yaml`)
- [x] Update `environment.yaml` with any new dependencies (scipy already available; netCDF4, eccodes already in env) (scipy, xarray, cfgrib if needed)

---

## Execution Order

```
1. OCP-Tool existing steps 0–3 (grid, FESOM mesh, LSM adaptation)
2. paleo_input.py: masks → ice → lakes → soils → soil water → vegetation → remaining fields → merge
3. paleo_topo.py: topography anomaly → spectral → write ICMSH
4. paleo_subgrid_oro.py: prepare → calnoro → post-process → apply to ICMGG
5. OCP-Tool existing steps 4+ (OASIS, runoff, etc.)
```

## Notes

- All interpolation done in Python (scipy) — no CDO dependency
- calnoro kept as external Fortran binary (called via subprocess)
- Reconstruction files follow PlioMIP3/PRISM4 conventions
- Lookup tables are hardcoded but configurable via config if needed later
