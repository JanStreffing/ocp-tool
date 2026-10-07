"""
OCP-Tool: OpenIFS Coupling Preparation Tool

Modules:
- config: Configuration loading and dataclasses
- cycles: OpenIFS cycle definitions (snow-field layout per cycle)
- grids: Gaussian grid processing
- lsm: Land-sea mask processing
- oasis_writer: OASIS file generation
- runoff: Runoff map modifications
- plotting: Visualization
- co2_interpolation: 3D CO2 interpolation
- field_interpolation: 2D field interpolation
- paleo_input: Paleo land surface modifications (ice, lakes, soils, vegetation)
- paleo_topo: Paleo topography modification (anomaly method)
- paleo_subgrid_oro: Subgrid-scale orography via calnoro
"""

# pyproj has to be loaded before ecCodes. The ecCodes wheel (eccodeslib) and the
# pyproj wheel each bundle native libraries, and when ecCodes is loaded first
# their teardown conflicts: Python aborts at exit with "double free or
# corruption" (exit code 134), after all the work is done. cartopy pulls in
# pyproj and lsm pulls in ecCodes, so without this the outcome depends on which
# ocp_tool module happens to be imported first.
try:
    import pyproj  # noqa: F401
except ImportError:
    pass

__version__ = "2.1.0"
