
__version__='0.1.1'

# raise custom exceptions onto the pytram package level
from .estimator import ExpressionError, NotConvergedWarning

# raise the dTRAM estimator class onto the pytram package level
from .dtram import DTRAM

# raise the API function onto the pytram level
from .api import dtram
