import os
from .CEA import *
from .utils import *
from .thermo_lib import printSimpleThermoLibLine

with open(os.path.join(os.path.dirname(__file__), "version.txt"), "r") as f:
    __version__ = f.read().strip()