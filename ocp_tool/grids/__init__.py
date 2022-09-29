from .regular import RegularLatLonGrid, FullGaussianGrid
from .orca import ORCA
from .gaussian import ReducedGaussianGrid
from .oifs import F128, TL159, TCO95, TL255, TCO159, TCO199, TCO319


def factory(grid_name, *args, **kwargs):

    reduced_gaussian_grids = {
        'TL159': TL159,
        'TCO95': TCO95,
        'TL255': TL255,
        'TCO159': TCO159,
        'TCO199': TCO199,
        'TCO319': TCO319,
    }

    full_gaussian_grids = {
        'F128': F128,
    }

    orca_grids = {
        'orca': ORCA,
        'ORCA': ORCA,
    }

    regular_latlon_grids = {
        'regular_latlon': RegularLatLonGrid,
    }

    if grid_name in reduced_gaussian_grids:
        return ReducedGaussianGrid(
            lats=reduced_gaussian_grids[grid_name].yvals,
            nlons=reduced_gaussian_grids[grid_name].reducedpoints,
        )

    elif grid_name in full_gaussian_grids:
        return FullGaussianGrid(
            # Note that we want a full Gaussian grid with latitudes starting at
            # the South pole, hence we have to reverse the FXXX lats
            lats=full_gaussian_grids[grid_name].yvals[::-1],
        )

    elif grid_name in orca_grids:
        return ORCA(*args, **kwargs)

    elif grid_name in regular_latlon_grids:
        return RegularLatLonGrid(*args, **kwargs)

    raise NotImplementedError(f'Unknown grid type: {grid_name}')
