"""
OpenIFS model-cycle definitions for OCP-Tool.

The layout of the surface fields in the ICMGG initial file changed between
OpenIFS cycles.  Everything that depends on the cycle is collected here so
that the rest of the code can stay cycle-agnostic.

The difference that matters for the paleo modifications is the snow scheme:

CY43R3 and earlier — single-layer snow (``LESNML = .false.``), GRIB1, surface::

     141  sd     snow depth              m of water equivalent
      32  asn    snow albedo             0-1
      33  rsn    snow density            kg m-3
     238  tsn    snow temperature        K

CY48R1 — multi-layer snow (``LESNML = .true.``, 5 layers), GRIB2,
``typeOfLevel = snowLayer``::

  228141  sd     snow mass               kg m-2
      32  asn    snow albedo             0-1
      33  rsn    snow density            kg m-3
     238  tsn    snow temperature        K
  228038  lwcs   snow liquid water       kg m-2

Note the unit change of the snow mass field: an ice-sheet point carrying
10 m of water equivalent is ``10.0`` in the CY43R3 field and ``10000.0`` in
the CY48R1 field.  Getting this wrong is silent, which is why the scale
factor lives in the cycle definition instead of in the calling code.

A CY48R1 executable with ``LESNML = .true.`` aborts during surface-field
initialisation::

    GRIDPOINT 3D FIELD MISSING:  228141  1
    ABOR1 : IOSTREAM_MIX:GRID_IN - MISSING FIELD

when it is handed a CY43R3 ICMGG file, because the 5-layer snow block is
absent.  The remaining CY48R1 additions (129172 C4 photosynthesis map,
200199 urban cover, 117-120 bare-soil albedo) are not generated here: they
are either carried through unchanged from the input file or appended by the
standard pipeline steps, so using a CY48R1 file as input is sufficient.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import eccodes


# ``typeOfLevel`` marking a message as belonging to the multi-layer snowpack.
SNOW_LAYER_TYPE_OF_LEVEL = 'snowLayer'

# Snow water equivalent assigned to a reconstruction ice-sheet point, in
# metres of water equivalent.  Scaled to the file units via CycleSpec.
ICE_SHEET_SNOW_MWE = 10.0


@dataclass(frozen=True)
class CycleSpec:
    """Cycle-dependent properties of the OpenIFS ICMGG initial file."""

    name: str
    snow_multi_layer: bool
    # GRIB code of the snow mass / snow depth field
    snow_mass_code: int
    # Multiply metres of water equivalent by this to obtain the file units
    snow_mass_scale: float
    snow_mass_units: str
    # GRIB codes carried once per snow layer (empty for single-layer cycles)
    snow_layer_codes: Tuple[int, ...] = ()

    def scale_snow_mass(self, metres_water_equivalent: float) -> float:
        """Convert m of water equivalent to the snow mass units of this cycle."""
        return metres_water_equivalent * self.snow_mass_scale

    @property
    def ice_sheet_snow_mass(self) -> float:
        """Ice-sheet snow mass in the units used by this cycle's ICMGG file."""
        return self.scale_snow_mass(ICE_SHEET_SNOW_MWE)


CY43R3 = CycleSpec(
    name='43r3',
    snow_multi_layer=False,
    snow_mass_code=141,
    snow_mass_scale=1.0,
    snow_mass_units='m of water equivalent',
)

CY48R1 = CycleSpec(
    name='48r1',
    snow_multi_layer=True,
    snow_mass_code=228141,
    snow_mass_scale=1000.0,
    snow_mass_units='kg m-2',
    # sd, asn, rsn, tsn, lwcs
    snow_layer_codes=(228141, 32, 33, 238, 228038),
)

CYCLES: Dict[str, CycleSpec] = {
    CY43R3.name: CY43R3,
    CY48R1.name: CY48R1,
}

# Spellings accepted in the YAML config
_ALIASES = {
    'cy43r3': '43r3', '43r3': '43r3', '43': '43r3',
    'cy48r1': '48r1', '48r1': '48r1', '48': '48r1',
}

AUTO = 'auto'


def normalise_cycle_name(name: str) -> str:
    """Normalise a user-supplied cycle name; returns ``'auto'`` unchanged."""
    key = str(name).strip().lower()
    if key in (AUTO, ''):
        return AUTO
    if key in _ALIASES:
        return _ALIASES[key]
    raise ValueError(
        f"Unknown OpenIFS cycle '{name}'. "
        f"Expected one of {sorted(CYCLES)} or '{AUTO}'."
    )


def get_cycle(name: str) -> CycleSpec:
    """Look up a CycleSpec by name (aliases accepted)."""
    key = normalise_cycle_name(name)
    if key == AUTO:
        raise ValueError("get_cycle() needs an explicit cycle, not 'auto'")
    return CYCLES[key]


def has_snow_layers(grib_file: Union[str, Path]) -> bool:
    """True if the GRIB file contains multi-layer snow (``snowLayer``) messages."""
    with open(grib_file, 'rb') as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                return False
            try:
                if eccodes.codes_get(gid, 'typeOfLevel') == SNOW_LAYER_TYPE_OF_LEVEL:
                    return True
            except eccodes.KeyValueNotFoundError:
                pass
            finally:
                eccodes.codes_release(gid)


def detect_cycle(grib_file: Union[str, Path]) -> CycleSpec:
    """
    Determine the cycle of an ICMGG file from its snow fields.

    Multi-layer snow (``typeOfLevel = snowLayer``) is present from CY48R1 and
    absent in CY43R3, which is the distinction the paleo code cares about.
    """
    return CY48R1 if has_snow_layers(grib_file) else CY43R3


def resolve_cycle(
    requested: Optional[str],
    grib_file: Optional[Union[str, Path]] = None,
    verbose: bool = True,
) -> CycleSpec:
    """
    Resolve the cycle to use, cross-checking config against the input file.

    Parameters
    ----------
    requested : cycle name from the config, or ``'auto'`` / ``None``
    grib_file : input ICMGG file used for detection / verification
    verbose : print the outcome

    Raises
    ------
    ValueError
        If an explicitly requested cycle contradicts the input file.  This is
        the mistake that otherwise surfaces much later as
        ``IOSTREAM_MIX:GRID_IN - MISSING FIELD`` at model start-up, or worse,
        as a silently wrong snow field.
    """
    key = normalise_cycle_name(requested if requested is not None else AUTO)
    file_exists = grib_file is not None and Path(grib_file).exists()

    if key == AUTO:
        if not file_exists:
            raise ValueError(
                "atmosphere.model_cycle is 'auto' but the input ICMGG file "
                f"({grib_file}) does not exist, so the cycle cannot be "
                "detected. Set atmosphere.model_cycle explicitly."
            )
        cycle = detect_cycle(grib_file)
        if verbose:
            print(f"  OpenIFS cycle: {cycle.name} (auto-detected from {grib_file})")
        return cycle

    cycle = CYCLES[key]
    if file_exists:
        detected = detect_cycle(grib_file)
        if detected.name != cycle.name:
            raise ValueError(
                f"Config requests OpenIFS cycle {cycle.name}, but {grib_file} "
                f"looks like {detected.name}: multi-layer snow messages are "
                f"{'absent' if cycle.snow_multi_layer else 'present'}. "
                f"Point paleo.icmgg_input_file at a {cycle.name} ICMGG file, "
                f"or set atmosphere.model_cycle to {detected.name}."
            )
    if verbose:
        print(f"  OpenIFS cycle: {cycle.name} (from config)")
    return cycle
