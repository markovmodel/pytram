
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# raise custom exceptions onto the pytram package level
from .estimator import ExpressionError, NotConvergedWarning

# raise the TRAMData class onto the pytram package level
from .reader import Reader

# raise the TRAMData class onto the pytram package level
from .tramdata import TRAMData

# raise the dTRAM estimator class onto the pytram package level
from ._dtram import DTRAM

# raise the xTRAM estimator class onto the pytram package level
from ._xtram import XTRAM

# raise the API function onto the pytram level
from .api import dtram, dtram_from_matrix, xtram, xtram_from_matrix

import warnings
warnings.warn("pytram is no longer supported; we recommend to use the PyEMMA package instead.")
del warnings
