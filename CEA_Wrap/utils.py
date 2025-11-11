import hashlib
import importlib.resources
import logging
import os
import platform
import shutil
from collections.abc import MutableMapping
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any, Iterator
from zlib import crc32

from .utils_low import getenv_t_f

log = logging.getLogger(__name__)

local_assets_directory = os.getenv("CEA_ASSETS_DIR") # or None
# If possible, use local assets rather than site-packages. But let the user force us to use just site-packages
# Gets set to false if there is some error with using local site packages
use_site_packages = getenv_t_f("CEA_USE_SITE_PACKAGES", False)

try_move_to_local = True # Whether or not we should try moving things from site packages
if use_site_packages: # If they specify to use site packages, they probably don't want to create a local directory
  try_move_to_local = False
else:
  if local_assets_directory is None:
    import appdirs  # Only import this if they don't specify a directory specifically
    local_assets_directory = appdirs.user_data_dir(__package__, roaming=False)  # By default we use the user's app data directory
  else:
    try_move_to_local = False # If the user sets the environment variable, we assume they don't want to move anything from an installation directory
    log.info(f"Local assets directory specified by environment variable! Using specified local assets directory: '{local_assets_directory}'")
    if not os.path.isdir(local_assets_directory):
      raise ImportError("Error! local assets directory specified by environment variable does not exist. Cannot start")

def get_site_package_path(file: str) -> str:
  # The reason the manager is used is because our package may be zipped and the manager extracts it
  #   However, this package is not zip-safe so we just return the location
  with importlib.resources.path(__package__+".assets", file) as manager:
    return str(manager)
    
def get_local_data_file(file: str) -> str:
  # Returns a file location in our data directory
  return os.path.join(local_assets_directory, file)

def get_asset(file: str) -> str:
  # Revert to package assets is required. Otherwise, get local data file if available.
  if use_site_packages:
    return get_site_package_path(file)
  else:
    return get_local_data_file(file)

"""
# The first time we install from source, we need to move our files from the ".assets" directory to our data directory
# This checks if the assets directory still exists, and if so, will try to copy files to it.
Expected process:
1st run: There is no AppData folder, all assets are in site-packages/CEA_Wrap/assets folder
          We create AppData/Local/CEA_Wrap and move all our assets there.
Subsequent runs: We check that the site-packages/CEA_Wrap/assets folder exists and it doesn't so we don't change anything
Subsequent installs (with --update): The assets folder exists, so we move all assets without replacement,
        so user thermo_spg.inp files are saved, but any missing assets are added
"""
if try_move_to_local: # If using pre-installed assets, we don't want to copy anything.
  try:
    log.info("Performing asset cleanup")
    asset_dir = get_site_package_path("")
    data_dir = get_asset("")
    if os.path.isdir(asset_dir):
      log.debug("Doing assets existence moving: package assets directory exists, moving files to data directory")
      log.debug("Package Dir: " + asset_dir)
      log.debug("Data Dir:    " + data_dir)
      if not os.path.isdir(data_dir): # Create our destination directory if it doesn't exist
        log.debug("Data directory didn't exist")
        os.makedirs(data_dir)
      for file in os.listdir(asset_dir):
        if file == "__pycache__": # Evidently the assets thing creates a pycache when it looks for paths
          continue
        src_path = os.path.join(asset_dir, file)
        dst_path = os.path.join(data_dir, file)
        log.debug("Checking file: " + file)
        log.debug(src_path + " ==> " + dst_path)
        if not os.path.exists(dst_path): # If we don't already have a copy of this file
          log.debug("Copying file")
          shutil.copy2(src_path, dst_path) # Copy it
        elif Path(dst_path).name == "trans.lib": # Special case for updating the trans lib
          log.debug("Examining trans.lib for hash")
          # Check if the existing destination file is the 12-byte "blank" trans.lib
          #   We want to check for this one exactly, in case the user has created their own trans.lib, we don't want to overwrite it
          # This is the hash of the "blank" file you get from no trans.inp
          BLANK_FILE_HASH = bytes.fromhex("E731531F76E6293EA5611B2F3B971A67AB7B1420")
          with open(dst_path, "rb") as f:
            if hashlib.sha1(f.read(), usedforsecurity=False).digest() == BLANK_FILE_HASH:
              log.debug("Detected legacy trans.lib. Overwriting with new version")
              shutil.copy2(src_path, dst_path)
            else:
              log.debug("trans.lib does not match legacy file hash. Not copying.")
        else:
          log.debug("File already existed in destination")
  except ModuleNotFoundError: # importlib raises this if a directory doesn't exist
    log.debug("Assets directory doesn't exist")
  except PermissionError: # Unsure why this happens, happened to one person using Anaconda
    log.warning("Permission Error. Cannot create local assets directory. Custom thermo.inp library will be overwritten when CEA_Wrap is updated")
    use_local_assets = False

def get_CEA_location(legacy=False) -> str:
  """
  Returns the location of the CEA executable, which may be overridden by the "CEA_EXE_LOCATION" environment variable

  :param legacy: If legacy, finds "FCEA2", otherwise finds "Dan_CEA2", defaults to False
  :return: The location of CEA
  """
  base_name = "FCEA2" if legacy else "Dan_CEA2"
  return os.getenv("CEA_EXE_LOCATION", get_asset(f"{base_name}.exe" if platform.system() == "Windows" else base_name))

def get_lib_locations() -> tuple[str, str]:
  """
  Returns the locations of the thermo.lib and trans.lib, which may be overridden by the environment variables
    "CEA_THERMO_LIB" and "CEA_TRANS_LIB", respectively

  :return: 2-tuple of locations of thermo.lib and trans.lib
  """
  return (
    os.getenv("CEA_THERMO_LIB", get_asset("thermo.lib")),
    os.getenv("CEA_TRANS_LIB",  get_asset("trans.lib"))
  )

def get_thermo_inp_location() -> str:
  """
  Returns the location of the thermo.inp input file, which is used to generate thermo.lib.
  The location may be overridden by the environment variable "CEA_THERMO_INP"

  :return: 
  """
  return os.getenv("CEA_THERMO_INP", get_asset("thermo_spg.inp"))

class Output(dict):
  """ This is just a dictionary that you can also use dot notation to access """
  def __getattr__(self, name) -> Any:
    return self[name]
  
  def __setattr__(self, name, value):
    if name.startswith("_"):
      super().__setattr__(name, value)
    else:
      self[name] = value

@dataclass
class DictDataclass(MutableMapping):
    """A dataclass that also supports dictionary-like access and methods."""
    
    def __getitem__(self, key: str) -> Any:
        """Allow dictionary-style access: obj['field']"""
        if not hasattr(self, key):
            raise KeyError(key)
        return getattr(self, key)
    
    def __setitem__(self, key: str, value: Any) -> None:
        """Allow dictionary-style assignment: obj['field'] = value"""
        if not hasattr(self, key):
            raise KeyError(key)
        setattr(self, key, value)
    
    def __delitem__(self, key: str) -> None:
        """Allow dictionary-style deletion: del obj['field']"""
        if not hasattr(self, key):
            raise KeyError(key)
        delattr(self, key)
    
    def __iter__(self) -> Iterator[str]:
        """Enable iteration over fields"""
        return iter(f.name for f in fields(self))
    
    def __len__(self) -> int:
        """Return number of fields"""
        return len(fields(self))
    
    def items(self):
        """Return (key, value) pairs"""
        return ((f.name, getattr(self, f.name)) for f in fields(self))
    
    def keys(self):
        """Return field names"""
        return (f.name for f in fields(self))
    
    def values(self):
        """Return field values"""
        return (getattr(self, f.name) for f in fields(self))
    
    def to_dict(self) -> dict[str, Any]:
      return {k:v for k, v in self.items()}

def move_file_if_changed(source: str, destination: str):
  """
  Copies the file at `source` to `destination` if a) `destination` does not exist or b) the CRC32 of `destination` does not match that of `source`

  :param source: The file path to the source file
  :param destination: The file path to the destination file
  :return: the destination
  """
  # file is the local destination, pack_file is the master location
  if os.path.isfile(destination):
    # If the file is here, we check the hash of it against the package one
    with open(source, "rb") as f1, open(destination, "rb") as f2:
      pack_hash = crc32(f1.read())
      local_hash = crc32(f2.read())
    if pack_hash != local_hash: 
      log.info(destination+" hash does not match package file hash! Updating local file with one from package")
      shutil.copyfile(source, destination)
  else:
    # If not here, copy it from package
    log.info(destination+" not found in current directory. Copying from package to current directory...")
    shutil.copyfile(source, destination)
  return destination

def open_thermo_lib():
  """ Opens the thermo library input file using the user's default .inp file viewer (should prompt if none) """
  print("Opening thermo library .inp file using default .inp file viewer")
  os.system(get_asset("thermo_spg.inp"))

def open_pdfs():
  """ Opens the NASA pdfs using the user's default pdf viewer """
  print("Opening user manuals using default pdf viewer")
  os.system('"'+get_asset("CEA_Mathematical_Analysis.pdf")+'"')
  os.system('"'+get_asset("CEA_Users_Manual_and_Program_Description.pdf")+'"')
  
def print_assets_directory():
  """ Just prints the directory where resources are """
  var = os.path.dirname(get_CEA_location())
  print(var)
  return var

###############################################################
################### Data Helper Utilities #####################
###############################################################

def _mutually_exclusive(*args):
  join_set = set()
  for arg in args:
    # First if this contains any elements that were in previous
    if not join_set.isdisjoint(arg):
      return False
    # Then update the current list of existing keys
    join_set.update(arg)
  return True # If none are combined, pass

class DataCollector(Output):
  def __init__(self, *args, keys=[], chamber_keys=[], throat_keys=[], exit_keys=[]):
    """
    DataCollector objects are Output objects that accepts lists of keys to make into lists
    :param keys: Also accepts list of arguments, these are keys such as 'cond' or 't_cp' or 'c_p'
    :param chamber_keys: List of chamber species mol/mass fractions to be included in the output object. Ex: "H2O" or "CO2". The key in the output will be the molecule name with "c_" prepended
    :param throat_keys: List of nozzle throat species mol/mass fractions to be included. The key in the output will be the molecule name with "t_" prepended
    :param exit_keys: List of nozzle exit species mol/mass fractions to be included. The key in the output will be the molecule name with nothing prepended.
    """
    if len(args) > 0:
      if keys:
        raise TypeError("Data_Collector should not receive both a list of arguments and the 'keys' keyword")
      keys = list(args)
    self._independent_variable=list()
    self._keys = keys
    self._chamber_keys = chamber_keys
    self._throat_keys = throat_keys
    self._exit_keys = exit_keys
    self._internal_keys = keys + exit_keys + list(map(lambda x: "t_"+x, throat_keys)) + list(map(lambda x: "c_"+x, chamber_keys))
    
    if not _mutually_exclusive(chamber_keys, throat_keys, exit_keys):
      raise ValueError("Can't have product keys in multiple keys arrays")
      
    if not all([isinstance(val, str) for val in keys + chamber_keys + throat_keys + exit_keys]):
      raise ValueError("All Data_Collector keys must be strings")
      
    for key in self._internal_keys:
      self[key] = list()
  
  @property
  def independent_variable(self):
    return self._independent_variable.copy()
  
  def sort(self):
    if len(self.independent_variable) != len(self[self._keys[0]]):
      raise ValueError("Independent variable array is not the same size as data arrays")
     # Sort indices by values: Enumerate returns iterator of (index, value). We then sort those by value (ascending). Make a list of the sorted indices
    index_order = [index for index, value in sorted(enumerate(self.independent_variable), key=lambda item: item[1])]
    # Re-order each internal list by the sorted index order
    for key in self._internal_keys:
      # https://stackoverflow.com/a/3260459
      self[key] = [self[key][index] for index in index_order]

  def add_data(self, data: Output|dict[str,Any], independent_variable=None):
    def try_add(key, inner):
      try:
        self[key].append(data[inner][key])
      except KeyError:
        self[key].append(0)

    if independent_variable is not None:
      self._independent_variable.append(independent_variable)

    for key in self._keys:
      self[key].append(data[key])
    for key in self._chamber_keys: # First go through chamber products
      try_add("c_"+key, "prod_c")
    for key in self._throat_keys: # Then go through throat products
      try_add("t_"+key, "prod_t")
    for key in self._exit_keys: # Then go through exit products
      try_add(key, "prod_e")
      
  def to_csv(self, filename: str, keys:list[str]|None=None, formatString="f"):
    """
    Writes all data to a csv, with the specified keys. If keys not given, uses all keys in order
    :param filename: The path at which to open a csv and write to it. Overwrites existing files
    :param keys: If None, will use all keys in object. Otherwise, a sorted list.
    :param formatString: The string to python's string formatting utility for each value. Default 'f'. Could be like "10.4f" or similar
    """
    log.info("Writing to CSV!")
    if keys is None:
      keys = sorted(self._internal_keys)
    
    with open(filename, "w") as file:
      if len(keys) == 0:
        log.warning("Tried writing DataCollector object to CSV with no keys")
        return
      for key in keys:
        file.write(key+",")
      file.write("\n")
      length = len(self[keys[0]])
      if length == 0:
        log.debug("Writing DataCollector object to CSV with no values")
      for i in range(length):
        for key in keys:
          file.write(("{:"+formatString+"},").format(self[key][i]))
        file.write("\n")

try:
  import numpy as np
  class NumpyDataCollector(DataCollector):
    def __init__(self, shape, *args, **kwargs):
      if not isinstance(shape, (list, tuple)):
        raise ValueError("shape parameter to NumpyDataCollector must be tuple of output shape")
      self._shape = shape
      
      super().__init__(*args, **kwargs)
      
      # Then we replace all the empty lists with zero-initialized arrays
      for key in self:
        self[key] = np.zeros(shape)
        
    def add_data(self, index, data):
      # index must be either the index for a 1-D array or a tuple for multi-dimensional arrays
      def try_add(key, inner):
        try:
          self[key][index] = data[inner][key]
        except KeyError:
          self[key][index] = 0
    
      for key in self._keys:
        self[key][index] = data[key]
      for key in self._chamber_keys: # First go through chamber products
        try_add(key, "prod_c")
      for key in self._throat_keys: # Then go through throat products
        try_add(key, "prod_t")
      for key in self._exit_keys: # Then go through exit products
        try_add(key, "prod_e")
except ImportError:
  import warnings
  warnings.warn("Numpy is not installed, 'NumpyDataCollector' module will not be available")
