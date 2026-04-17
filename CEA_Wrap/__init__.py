import os as _os
from .wrapper import *
from .utils import *
from .thermo_lib import *

with open(_os.path.join(_os.path.dirname(__file__), "version.txt"), "r") as f:
    __version__ = f.read().strip()
