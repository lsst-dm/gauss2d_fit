import pkgutil

__path__ = pkgutil.extend_path(__path__, __name__)

import lsst.gauss2d

from ._gauss2d_fit import *
from ._gauss2d_fit import __version__
