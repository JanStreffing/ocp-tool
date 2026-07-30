"""OASIS remapping weight generation (thin compatibility wrapper).

The real implementation now lives in :mod:`ocp_tool.oasis_weights`, which drives
the coupled model's own pyOASIS (one MPI rank per coupling link) to compute the
``rmp_*.nc`` weight files. This module keeps the ``generate_rmp_files`` entry
point that ``run_ocp_tool.py`` calls, delegating to the new engine.
"""

from pathlib import Path
from typing import Dict, List, Optional

from ocp_tool.oasis_weights import Link, awiesm3_feom_links, generate_weights


def generate_rmp_files(
    oasis_dir,
    coupling_links: Optional[List[Dict[str, str]]] = None,
    component_name: str = "ocp_tool",
    *,
    method: str = "existing",
    atm_grid: str = "A096",
) -> List[Path]:
    """Generate OASIS remapping weight files in ``oasis_dir``.

    ``oasis_dir`` must already contain grids.nc / masks.nc / areas.nc with every
    grid referenced by the links (ocp-tool writes these; the feom grid for the
    current FESOM mesh must be present).

    ``coupling_links`` may be a list of dicts ({"source","target","scripr"
    [,"loctrans"]}); if omitted, the default awiesm3 feom links are used.
    ``method`` selects the remap profile ("existing" reproduces the current
    non-conservative namcouple; "conserv" uses feom corners for conservative
    flux remapping).
    """
    links = None
    if coupling_links:
        links = [
            lk if isinstance(lk, Link) else Link.from_dict(lk)
            for lk in coupling_links
        ]
    return generate_weights(
        oasis_dir,
        links=links,
        atm_grid=atm_grid,
        method=method,
    )


__all__ = ["generate_rmp_files", "awiesm3_feom_links", "Link", "generate_weights"]
