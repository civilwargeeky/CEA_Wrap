import os
import tempfile
import shutil
import subprocess
from typing import Optional, TextIO
import logging


from .utils import get_CEA_location, get_lib_locations, move_file_if_changed
from .utils_low import getenv_t_f

log = logging.getLogger(__name__)

DEFAULT_COPY_TO_TEMP_DIR = getenv_t_f("CEA_COPY_THERMO_TO_TEMP", True)

class CEA:
  """
  Objects of this class represent the locations of all files needed to run CEA and code for interacting with the CEA executable/backend
  """
  OUT_ERROR_FILE = "CEA_Wrap_error.out"
  DELIMITER = '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+' # Delimiter in CEA output between .out and .plt files
  CEA_LEGACY = False # Parameter given to get_CEA_location for getting CEA exe

  def __init__(self, CEA_path:Optional[str]=None, thermo_lib_path:Optional[str]=None, trans_lib_path:Optional[str]=None,
               copy_to_temp_dir:bool=DEFAULT_COPY_TO_TEMP_DIR, dump_out_on_error=True):
    """
    Create a new CEA interface object. If paths are not given, they are taken from the default assets directory.

    :param CEA_path: Path to the CEA executable. This executable should accept a .inp file from STDIN
                     and output both the .out and .plt files (separated by a delimiter) to STDOUT (not typical FCEA2 behavior), defaults to None
    :param thermo_lib_path: Path the thermo.lib file, defaults to None
    :param trans_lib_path: Path to the trans.lib file, defaults to None
    :param copy_to_temp_dir: If True, thermo.lib and trans.lib will be copied to a temporary directory on initialization. This is useful
                             when using CEA_Wrap in a multithreaded or multiprocessing context, as CEA.exe locks both files while running, defaults to True
    :param dump_out_on_error: If True, when CEA errors it will dump the .out file to "CEA_Wrap_error.out"
    :NOTE: If you are copying a temp dir, it's a good idea to treat this as a singleton so we don't have lots of copies of thermo.lib floating around.
    """
    
    self.CEA_path = get_CEA_location(self.CEA_LEGACY) if CEA_path is None else CEA_path
    default_thermo_lib, default_trans_lib = get_lib_locations()
    # We have both an "original" and "path" because we may move this to a temporary directory
    self.original_thermo_lib = self.thermo_lib_path = default_thermo_lib if thermo_lib_path is None else thermo_lib_path
    self.original_trans_lib  = self.trans_lib_path  = default_trans_lib  if trans_lib_path  is None else trans_lib_path

    self._copy_to_temp_dir = None # Local variable. Set to neither True nor False so cleanup gets handled immediately
    self.temporary_directory = None
    self.copy_to_temp_dir = copy_to_temp_dir # getter/setter that handles file creation

    self.dump_out_on_error = dump_out_on_error
  
  @property
  def copy_to_temp_dir(self) -> bool:
    return self._copy_to_temp_dir
  
  def _populate_temp_dir(self):
    # Note: This directory will be removed when this object is deleted, which happens when our CEA object is deleted
    log.debug("Creating temporary directory for thermo.lib files")
    self.temporary_directory = tempfile.TemporaryDirectory(prefix="CEA_Wrap_")
    log.debug(f"Temporary folder created at '{self.temporary_directory.name}'")
    self.thermo_lib_path = os.path.join(self.temporary_directory.name, "thermo.lib")
    self.trans_lib_path = os.path.join(self.temporary_directory.name, "trans.lib")
    for src, dst in ((self.original_thermo_lib, self.thermo_lib_path), (self.original_trans_lib, self.trans_lib_path)):
      log.debug(f"  Copying '{src}' to '{dst}'")
      shutil.copy(src, dst)
  
  def _cleanup_temp_dir(self):
    if self.temporary_directory is not None:
      log.debug(f"Cleaning up existing temporary directory in {self.temporary_directory.name}")
      self.temporary_directory.cleanup() # explicitly delete files and folder before removing reference (probably not needed)
      self.temporary_directory = None
    self.thermo_lib_path = self.original_thermo_lib # Reset our references
    self.trans_lib_path = self.original_trans_lib
  
  @copy_to_temp_dir.setter
  def copy_to_temp_dir(self, new_value: bool):
    """
    Changes the value of copy_to_temp_dir, creating a temp directory or deleting it if necessary
    """
    if new_value == self._copy_to_temp_dir:
      return
    if new_value: # We need to create a new temporary directory
      self._populate_temp_dir()
    else:
      self._cleanup_temp_dir()

  def run_cea_backend(self, *, contents:None|str=None, file:None|TextIO=None, filename:None|str=None) -> tuple[str, str]:
    """
    Actually runs the CEA executable with the .inp file contents, which may be specified as a string, file-like object, or file location (but only one of those)

    :param contents: The contents of a .inp file to pass to CEA
    :param file: A file-like object containing a .inp file
    :param filename: The path to a .inp file
    :raises RuntimeError: If CEA errors, raises RuntimeError
    :return: Returns the contents of the .out file and the contents of the .plt file
    """
    if filename is not None:
      with open(filename) as f:
        contents = f.read()
    elif file is not None:
      contents = file.read()
    contents = contents.replace("\n", "\r\n")

    # Path to program, with arguments of thermo.lib location and trans.lib location
    arguments = [self.CEA_path, self.thermo_lib_path, self.trans_lib_path]
    # Run the program, with file contents as input to stdin and capturing stdout
    # NOTE: We cannot pass the contents in as a text buffer via PIPE, because CEA seeks forward *and backward* in the file
    #       Default STDIN with a PIPE cannot seek backward, and results in undefined behavior.
    # The temporary file is deleted after CEA is run
    with tempfile.TemporaryFile(prefix="CEA_Wrap_") as io_file:
      io_file.write(contents.encode("utf-8"))
      io_file.seek(0)
      process = subprocess.run(arguments, stdin=io_file, text=True, capture_output=True)
    if process.returncode != 0:
      log.error(process)
      if self.dump_out_on_error:
        try:
          with open(self.OUT_ERROR_FILE, "w") as f:
            f.write(process.stdout)
        except Exception:
          raise RuntimeError(f"Running CEA failed with errors. Could not write output dump file")
        else: # If we write the output error dump file successfully
          raise RuntimeError(f"Running CEA failed with errors. Output dumped to {self.OUT_ERROR_FILE}")
      else:
        raise RuntimeError("Running CEA failed with errors")
        
    output = process.stdout # Contains both the .out and .plt files
    if self.DELIMITER not in output:
      raise ValueError("Delimiter not found in CEA output. Cannot process data. Are you using the correct executable?")
    out, plt = map(str.strip, output.split(self.DELIMITER)) # Split the output into two files, separated by the delimiter

    return out, plt

class Legacy_CEA(CEA):
  """ 
  Uses legacy file-writing and reading behavior with standard FCEA2 executable
  """
  CEA_LEGACY = True # Will get legacy CEA executable by default
  INFILE_PREFIX = "_CEA_calculation" # PREFIX for .inp, .out, .plt files used by CEA

  def _cleanup_temp_dir(self):
    """
    We can still use the temporary directory with trans.lib and thermo.lib created in the CEA class, but if they don't want a temp directory, we
    need to put everything in the current directory
    """
    super()._cleanup_temp_dir() # Still cleanup the temp directory if it exists
    # Then we need to move files into current directory so we can use them
    log.debug("Moving thermo lib and trans lib into current working directory")
    try:
      self.thermo_lib_path = move_file_if_changed(self.original_thermo_lib, os.path.join(os.getcwd(), "thermo.lib"))
      self.trans_lib_path  = move_file_if_changed(self.original_trans_lib , os.path.join(os.getcwd(), "trans.lib" ))
    except PermissionError as e:
      log.error("---- Error! Attempted to copy thermo.lib and trans.lib into current directory but failed ----")
      log.error("---- Is your current directory system32 or another protected directory? ----")
      raise e from None
  
  def run_cea_backend(self, *, contents:None|str=None, file:None|TextIO=None, filename:None|str=None) -> tuple[str, str]:
    """
    Actually runs the CEA executable with the .inp file contents, which may be specified as a string, file-like object, or file location (but only one of those)

    :param contents: The contents of a .inp file to pass to CEA
    :param file: A file-like object containing a .inp file
    :param filename: The path to a .inp file
    :raises RuntimeError: If CEA errors, raises RuntimeError
    :return: Returns the contents of the .out file and the contents of the .plt file
    """
    # Directory of all files to process
    directory = os.getcwd() if self.temporary_directory is None else self.temporary_directory.name
    # Expected location of all CEA files
    inp_file = os.path.join(directory, self.INFILE_PREFIX+".inp")
    out_file = os.path.join(directory, self.INFILE_PREFIX+".out")
    plt_file = os.path.join(directory, self.INFILE_PREFIX+".plt")
    
    if contents is not None:
      with open(inp_file, "w") as f:
        f.write(contents)
    elif file is not None:
      with open(inp_file, "w") as f:
        shutil.copyfileobj(file, f)
    elif filename is not None:
      shutil.copyfile(filename, inp_file)
    
    # Path to program, with arguments of thermo.lib location and trans.lib location
    arguments = [self.CEA_path, self.thermo_lib_path, self.trans_lib_path]
    # Run the program, executing it within the directory with thermo.lib and the .inp file
    # Pass "\n" after the prefix of the file so that the line is run
    process = subprocess.run([self.CEA_path], cwd=directory, input=self.INFILE_PREFIX+"\n", text=True, stdout=subprocess.DEVNULL)
    if process.returncode != 0:
      log.error(process)
      if self.temporary_directory is None:
        if self.dump_out_on_error:
          try:
            shutil.copyfile(out_file, self.OUT_ERROR_FILE)
          except Exception:
            raise RuntimeError(f"Running CEA failed with errors. Could not write output dump file")
          else: # If we write the output error dump file successfully
            raise RuntimeError(f"Running CEA failed with errors. Output dumped to {self.OUT_ERROR_FILE}")
        else: # If no dumping on errors
          raise RuntimeError("Running CEA failed with errors")
      else: # If we are just running in the current working directory
        raise RuntimeError(f"Running CEA failed with errors. Check output within {out_file}")
    
    try:
      with open(out_file) as f:
        out_text = f.read()
    except FileNotFoundError: # CEA always generates a .out file
      log.error(f".out file was not found in expected location: {out_file}")
      raise RuntimeError(".out file was not generated from run. Cannot process CEA output")
    
    try:
      with open(plt_file) as f:
        plt_text = f.read()
    except FileNotFoundError: # It is not required for a problem to use the .plt file, so it may not be generated
      plt_text = ""
    
    return out_text, plt_text