################################################################################
#
#   test_xtram.py - testing the xTRAM estimator
#
#   author: Antonia Mey <antonia.mey@fu-berlin.de>
#
################################################################################

from nose.tools import assert_raises, assert_true
from pytram import XTRAM, ExpressionError, NotConvergedWarning
import numpy as np


