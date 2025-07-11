"""
This module contains warning filters to silence specific warnings in the ocp-tool.
Import this module at the beginning of your script to suppress the warnings.
"""

import warnings

# Silence the Dask/distributed large graph size warnings
warnings.filterwarnings("ignore", message="Sending large graph of size .* MiB", 
                      category=UserWarning, module="distributed.client")
