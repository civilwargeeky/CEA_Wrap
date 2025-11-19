import logging
import os
import re
import threading
from dataclasses import dataclass, field
from typing import Callable, Generic, Never, Sequence, TypeVar, get_args, Literal

from .utils_low import getenv_t_f
from .utils import DictDataclass, Output

if getenv_t_f("CEA_USE_LEGACY", False):
  from .cea_interface import Legacy_CEA as CEA_Class
else:
  from .cea_interface import CEA as CEA_Class

from .thermo_lib import ThermoInterface as ThermoInterface_Class

log = logging.getLogger(__name__)

TempUnitType = Literal["k", "r", "c", "f"] # options for temperature units
PresUnitType = Literal["bar", "atm", "psi", "mmh"] # options for pressure units
EnthalpyUnitType = Literal["c", "kc", "j", "kj"] # options for enthalpy units. 'c' and 'kc' are "calorie" and "kilo-calorie (Calorie)"
OPTIONS_TEMP_UNITS = get_args(TempUnitType)
OPTIONS_PRES_UNITS = get_args(PresUnitType)
OPTIONS_ENTHALPY_UNITS = get_args(EnthalpyUnitType)

# Instantiate default objects used for all problems (using files defined within .utils.py)
# Because the CEA executable locks the thermo.lib and trans.lib files while running, only one instance of the executable can access those files at a time
#   to allow for multiprocessing, each instance of CEA_Class will copy those files to a temporary directory
#   To allow for threading, we could either a) use a global mutex on that subprocess call or b) create a temp directory per-thread
#   To increase performance, we will create a temp-directory per-thread by default
CEA_thread_dict = {
  threading.get_ident(): CEA_Class()
}
ThermoInterface = ThermoInterface_Class()

@dataclass
class ChemicalRepresentation:
  """
  A ChemicalRepresentation represents a custom material added in the .inp file, as opposed to using an existing material in thermo.lib

  :param chemical_representation: chemical composition in space-separated pairs of (element, n_atoms), such as "LI 1 B 1 H 4" for LiBH4.
  :param hf: Enthalpy of formation, in `hf_unit`/mol.
  :param hf_unit: Unit for enthalpy of formation. defaults to "kj"
  :param is_internal_energy: If True, sets an internal energy (`u`), rather than enthalpy (`h`). defaults to False
  """
  chemical_composition: str
  hf: float
  hf_unit: EnthalpyUnitType = "kj"
  is_internal_energy:bool = False
  
  def __post_init__(self):
    if self.hf_unit not in OPTIONS_ENTHALPY_UNITS:
      raise ValueError(f"hf_unit was '{self.hf_unit}' but must be one of {OPTIONS_ENTHALPY_UNITS}")
  
  def get_CEA_string(self):
    return f"{self.chemical_composition} {'u' if self.is_internal_energy else 'h'},{self.hf_unit}/mol={self.hf:.5f}"

# Function to let people know that something is deprecated
def deprecated_setter(fcn: Callable):
  def subfcn(*args, **kwargs):
    import warnings
    warnings.warn(f"The function {fcn.__name__} is deprecated. Please use the setter instead", DeprecationWarning)
    fcn(*args, **kwargs)
  subfcn.__name__ = fcn.__name__
  return subfcn


class Material:
  output_type = None # MUST BE SPECIFIED IN SUBCLASSES. The string we put before the material when writing .inp files
  is_fuel = None # MUST BE SPECIFIED IN SUBCLASSES. Assume all materials are either fuels or oxidizers. Must be overridden
  # If true, on initialization, we check if the material exists in the supplied thermo_spg input file
  #   also, will add .ref member to all objects for the thermo input reference
  check_against_thermo_inp = True
  
  def __init__(self, name, temp:float=298.15, wt:float|None=None, mols:float|None=None, chemical_representation:ChemicalRepresentation|None=None,
               wt_percent:Never=None):
    """
    A material to put in the problem's .inp file. Can be specified inline or come from the thermo.lib
    If neither wt nor mols is specified, wt is 100

    :param name: Name to input into CEA
    :param temp: temperature in K, defaults to 298.15K
    :param wt: Relative weight percent of this Material, defaults to None.
    :param mols: Mole ratio for this material, defaults to None
    :param chemical_representation: chemical composition, including name and hf. Such as `ChemicalRepresentation("LI 1 B 1 H 4", -200)` for LiBH4. 
        If you define a representation for a chemical already in the thermo library, this will ignore thermo library properties. Parameter defaults to None
    :raises ValueError: _description_
    :raises TypeError: _description_
    """
    if wt_percent is not None:
      import warnings
      warnings.warn("The parameter has been renamed 'wt'. Setting wt_percent is no longer supported and will be removed in a future release.", DeprecationWarning)
      wt = wt_percent
    if wt is None and mols is None: # If neither is specified, user probably doesn't care
      wt = 100
    elif wt is not None and mols is not None:
      raise TypeError("Material cannot have both wt and mols specified")
    
    self._name = name
    self._temp = temp
    self._wt = wt
    self._mols = mols
    self._chemical_representation = chemical_representation
    
    self.ref = None
    if self.check_against_thermo_inp and not chemical_representation: # If they don't override thermo.lib data and we check for missing elements
      if name not in ThermoInterface:
        close_matches = ThermoInterface.get_close_matches(name)
        raise ValueError(f"specified element '{name}' does not exist in package thermo library\n" +
                         f"Change name or set {__package__}.Material.check_against_thermo_inp to False\n"+
                         f"Material '{name}' {len(close_matches)} closest matches: \"" + '", "'.join(close_matches)+'"')
      
      self.ref = ThermoInterface[name] # Get the ThermoMaterial and store it in ref
    
      if not ThermoInterface[name].defined_at(temp):
        ranges = ThermoInterface[name].temp_ranges
        # min and max sort by 1st element of tuples, assume tuple with max 1st element is max total.
        string = f"specified material '{name}' does not exist at temperature {temp:0.2f} " +\
                  " (min, max) = ({:.2f}, {:0.2f})".format(min(ranges)[0], max(ranges)[1]) 
        raise ValueError(string)

  @property
  def name(self): # Name cannot be set outside of constructor
    return self._name
  
  @property
  def wt(self) -> float:
    if self._wt is None:
      raise ValueError("Material had no weight")
    return self._wt
  
  @wt.setter
  def wt(self, value:float):
    self._mols = None # Can only have one
    self._wt = value

  @property
  def mols(self) -> float:
    if self._mols is None:
      raise ValueError("Material had no mols")
    return self._mols
  
  @mols.setter
  def mols(self, value:float):
    self._wt = None # Can only have one
    self._mols = value

  @property
  def temp(self):
    return self._temp

  @temp.setter
  def temp(self, value:float):
    if self.check_against_thermo_inp and not ThermoInterface[self.name].defined_at(value):
      raise ValueError(f"specified material '{self.name}' does not exist at temperature {value:0.2f}")
    self._temp=value

  temperature = temp # Alias for "temp"

  ### set functions to explicitly set values for the Material

  def set_wt(self, value: float):
    self.wt = value # Call new setter

  def set_mols(self, value:float):
    self.mols = value # Call new setter
  
  def set_temp(self, value:float):
    self.temp = value # Call new setter
  
  def set_temperature(self, value: float):
    self.temp = value # Call underlying setter
  
  @property
  @deprecated_setter
  def chemical_composition(self):
    return self._chemical_representation # Cannot be set outside of constructor

  @property
  def chemical_representation(self):
    return self._chemical_representation # Cannot be set outside of constructor
    
  def is_mols(self) -> bool: # Helper function for a material being in wt or mols
    return self._mols is not None
    
  def is_wt(self) -> bool: # Helper function for a material being in wt or mols
    return self._wt is not None

  def get_CEA_name_wt(self, wt_or_mol: str, amount: float) -> str:
    """ Gives the actual line put into .inp files, minus potential hf strings or chemical composition. Not valid if self.composition is given """
    return f"   {self.output_type}={self.name} {wt_or_mol}={amount:0.5f}  t,k={self.temp:0.5f}"

  def get_CEA_str(self) -> str:
    # Specify whether to use the str/val for weight or mols.
    name = "wt" if self.is_wt() else "mol"
    ratio = self.wt if self.is_wt() else self.mols

    if ratio > 0:
      string = self.get_CEA_name_wt(name, ratio)
      if self._chemical_representation:
        string += self._chemical_representation.get_CEA_string()
      return string + "\n"
    elif ratio == 0: # This allows us to not include materials that are set to 0
      return ""
    else: # If ratio < 0
      raise ValueError("Cannot have {} with <0 {} percent!".format(self.name, name))
      
  def __str__(self):
    return self.name

  ### These functions for setting values are deprecated and will be removed in a future release
  @property
  @deprecated_setter
  def wt_percent(self):
    return self.wt

  @deprecated_setter
  def set_wt_percent(self, value:float):
    self.wt = value # Call new setter
      
class Fuel(Material):
  output_type = "fuel"
  is_fuel = True

F = Fuel # Alias

class Oxidizer(Material):
  output_type = "ox"
  is_fuel = False

O = Oxidizer # Alias

def make_fuel_property(name: str) -> property:
  """ Helper for Problem objects to make properties for all the different fuel ratio types """
  def getter(self: "Problem"): 
    f"""
    Get/set the "{name}" fuel property for this Problem
    If the current fuel property is not "{name}" and you try to get the value, a ValueError will be raised
    """
    if self.ratio_name == name:
      return self.ratio_value
    raise ValueError(f'Fuel ratio was "{self.ratio_name}", not {name}')
  
  def setter(self: "Problem", value: float):
    return self._set_fuel_ratio(**{name: value})
  
  return property(getter, setter)

OutputType = TypeVar("OutputType")

class Problem(Generic[OutputType]):
  ### NOTE : ALL PROBLEM SUBCLASSES MUST SPECIFY PROBLEM TYPE AND PLOT KEYS ###
  problem_type = None
  plt_keys: str # plt_keys must be defined by subclasses. It should be a space-separated string of items to go into the "plt" command of FCEA2
  
  # If true, on initialization, we check if all inserts and omits exist in the supplied thermo_spg input file
  check_against_thermo_inp = True 
  
  @classmethod
  def _get_plt_keys(cls) -> str:
    try:
      return cls.plt_keys
    except AttributeError:
      raise NotImplementedError(f"Class {cls.__name__} should implement .plt_keys, but did not") from None

  _ratio_options = ["p_f", "f_o", "o_f", "phi", "r_eq"] # Possible function arguments
  _ratio_CEA =     ["%f",  "f/o", "o/f", "phi", "r"] # Values to put into CEA
  def _set_fuel_ratio(self, **kwargs):
    # This function will go through kwargs and find which ratio we should use
    found_one = False
    for kwarg_name, cea_name in zip(self._ratio_options, self._ratio_CEA):
      if kwarg_name in kwargs:
        if found_one: # If we already found one, we can't use another
          raise TypeError("Can only specify one ratio at once (%f, o/f, phi, etc.)")
        self.ratio_name = cea_name # put in the string to place into CEA
        self.ratio_value = kwargs[kwarg_name] # Then get the float
        found_one = True
  
  def _format_input_list(self, input_list: str | list[str] | list[Material] | None):
    """
    Transmutes the input_list to a list of names
       checks each element to ensure it is a valid material

    :param input_list: Either space-separated str of material strings or a list of the same
    :raises ValueError: If we check_against_thermo_inp and a material isn't found in ThermoInterface, raises an error
    :return: List of validated material names
    """
    
    if input_list is None:
      return None
    if isinstance(input_list, str):
      input_list = input_list.split()
    input_list = list(map(str, input_list)) # make any Materials supplied into strings
    if self.check_against_thermo_inp:
      for material in input_list:
        if material not in ThermoInterface:
          close_matches = ThermoInterface.get_close_matches(material)
          raise ValueError(f"specified element '{material}' does not exist in package thermo library\n" +
                           f"Change name or set {__package__}.{__class__}.check_against_thermo_inp to False\n"+
                           f"Material '{material}' {len(close_matches)} closest matches: \"" + '", "'.join(close_matches)+'"')
    return input_list
    
  # All arguments must be specified by keyword
  def __init__(self, *, 
               pressure: float=1000, 
               materials: list[Material]|None=None, 
               massf: bool=False, 
               pressure_units: str="psi", 
               inserts: str | list[str] | list[Material] | None=None, 
               omits:   str | list[str] | list[Material] | None=None,
               filename:str|None=None,
               **kwargs
               ):
    """
    A generic CEA problem

    :param pressure: Chamber pressure, defaults to 1000
    :param materials: List of materials to run in this problem, defaults to None
    :param massf: Output product species as mass fractions, rather than mole fractions, defaults to False
    :param pressure_units: Units for given pressure, defaults to "psi"
    :param inserts: List of species to force into consideration for product species. Useful for problems with a lot of condensed phase in the products. 
                    Can be specified as a space-separated list of names, list of names, or list of materials, defaults to None
    :param omits: List of species to force out of consideration for product species. Can be specified as a space-separated list of names,
                  list of names, or list of materials, defaults to None, defaults to None
    :param [fuel_ratio]: Fuel ratio may be specified o_f, phi, or others (see full documentation)
    :param filename: Deprecated: No longer needed as the fixed CEA executable does not write to or read from files. Used to be the 
                     prefix for the .inp files to use for this problem, .out and .plt files will share this
    :raises TypeError: If extra parameters are set that are not one of the fuel_ratios
    """
      
    self.massf = massf
    self.materials: list[Material] = materials or []
    self.pressure = pressure
    self._pressure_units = None
    self.set_pressure_units(pressure_units)
    self.filename:str # Specify before setting
    self.set_filename(filename)
    
    # Check against thermo file for these to prevent errors
    self._inserts = self._format_input_list(inserts)
    self._omits   = self._format_input_list(omits)
    
    self.ratio_name = None
    self.ratio_value = None
    self._set_fuel_ratio(**kwargs)
    
    if (curr_thread := threading.get_ident()) in CEA_thread_dict: # Get the CEA object local to our thread
      self.CEA = CEA_thread_dict[curr_thread]
    else: # If no thread-local object, instantiate a new one
      self.CEA = CEA_thread_dict[curr_thread] = CEA_Class()
    
    diff = set(kwargs).difference(set(self._ratio_options)) # Find if kwargs contains any entries not in self._ratio_options
    if diff: # If there are any kwargs keys that aren't in _ratio_options keys, we should error
      raise TypeError(self.__class__.__name__+"() got an unexpected keyword argument(s): " + ",".join(diff))
  
  # This sets up "properties" for all the fuel properties
  # This allows people to do `my_problem.f_o = 5` or `problem.phi = 1`
  p_f = make_fuel_property("p_f")
  f_o = make_fuel_property("f_o")
  o_f = make_fuel_property("o_f")
  phi = make_fuel_property("phi")
  r_eq = make_fuel_property("r_eq")

  # Setter functions for all fuel ratios
  def set_p_f(self, p_f: float): self._set_fuel_ratio(p_f=p_f)
  def set_f_o(self, f_o: float): self._set_fuel_ratio(f_o=f_o)
  def set_o_f(self, o_f: float): self._set_fuel_ratio(o_f=o_f)
  def set_phi(self, phi: float): self._set_fuel_ratio(phi=phi)
  def set_r_eq(self, r_eq: float): self._set_fuel_ratio(r_eq=r_eq)

  # Set up properties for inserts and omits
  @property
  def inserts(self): return self._inserts
  @inserts.setter
  def _insert_setter(self, value: str | list[str] | list[Material] | None): return self._format_input_list(value)

  @property
  def omits(self): return self._omits
  @omits.setter
  def _omit_setter(self, value: str | list[str] | list[Material] | None): return self._format_input_list(value)
  
  def set_pressure(self, pressure: float): self.pressure = pressure
  def set_materials(self, materials: list[Material]): self.materials = materials
  def set_massf(self, massf: bool): self.massf = massf
  def set_inserts(self, inserts: str | list[str] | list[Material] | None): self._inserts = self._format_input_list(inserts)
  def set_omits(self, omits: str | list[str] | list[Material] | None): self._omits = self._format_input_list(omits)

  def set_filename(self, filename: str|None):
    if filename is not None:
      import warnings
      warnings.warn("setting filename is no longer needed for CEA problems and this key is ignored", DeprecationWarning)
  
  @property
  def pressure_units(self): return self._pressure_units
  @pressure_units.setter
  def _pressure_units_setter(self, value: str): return self.set_pressure_units(value)

  def set_pressure_units(self, pressure_units: str):
    pressure_units = pressure_units.lower() # We only have a few options for pressure units
    if pressure_units not in OPTIONS_PRES_UNITS:
      raise ValueError("pressure unit must be in " + ", ".join(OPTIONS_PRES_UNITS))
    self._pressure_units = pressure_units
  
  def set_absolute_o_f(self) -> float:
    # Set an o_f ratio assuming that the wt of each of our materials is an absolute percentage
    sum_ox   = sum([item.wt for item in filter(lambda x: not x.is_fuel, self.materials)], start=0.0)
    sum_fuel = sum([item.wt for item in filter(lambda x: x.is_fuel, self.materials)], start=0.0)
    o_f = sum_ox/sum_fuel
    self.set_o_f(o_f)
    return o_f
  
  def run(self, *materials:Material) -> OutputType:
    if self.ratio_name == None:
      raise TypeError("No reactant ratio specified, must set phi, or o/f, or %f, etc.")
    if len(materials) > 0: # If they specify materials, update our list
      self.materials = list(materials)
    
    try:
      file_contents = self.make_input_file(self.materials)
    except OSError: # Sometimes things like dropbox lock the file so we can't access it
      raise RuntimeError("unable to open input file for writing...")
    
    out_file, plt_file = self.CEA.run_cea_backend(contents=file_contents)
    return self.process_output(out_file, plt_file)
  
  def run_cea(self, *args, **kwargs): # alias for backward-compatibility
    import warnings
    warnings.warn("run_cea is deprecated. Use 'run' instead", DeprecationWarning)
    return self.run(*args, **kwargs)
  
  def make_input_file(self, material_list: Sequence[Material]) -> str: # chamber conditions and materials list
    # Make sure we have some materials
    if material_list is None or len(material_list) == 0:
      raise ValueError("must specify at least one Material in Problem")

    # Makes sure that all our materials are Material objects
    if any([not isinstance(x, Material) for x in material_list]):
      raise ValueError("all items in Problem material list must be Material objects")
    # First we need to check that all of our materials have either mole/wt fractions specified. No mixing.
    # If the first element has a wt, then we check that all elements have a wt
    if not all([(mat.is_mols() if material_list[0].is_mols() else mat.is_wt()) for mat in material_list]):
      raise ValueError("all materials must use wt percent or mol ratio, not a mixture of both")

    if any((type(x) == Material) for x in material_list):
      raise ValueError("tried running problem with base Material. All Material inputs must be Fuel or Oxidizer")

    fuels = list(filter(lambda x: x.is_fuel, material_list))
    oxidizers = list(filter(lambda x: not x.is_fuel, material_list))
    
    # Here we check for monopropellant cases. In a monopropellant case, you must specify only fuels and an o_f of 0 or you get weird results
    if len(fuels) == 0 and len(oxidizers) == 1:
      raise ValueError("Monopropellant problems must only specify fuels, not oxidizers")
    if len(oxidizers) == 0 and (self.ratio_name != "o/f" or self.ratio_value != 0):
      raise ValueError("Monopropellant/all-fuel problems must be run with o/f of 0 (set o_f=0)")
    if len(fuels) == 0:  # otherwise, if we have no fuels but multiple oxidizers, more general error
      raise ValueError("Must specify at least one fuel")
  
    contents = ""
    contents += "problem "
    contents += self.get_prefix_string()
    contents += "react  \n"
    # First we write all the fuels (I don't think there's a reason to write them in this order, but we'll do it anyway
    for fuel in fuels:
      contents += fuel.get_CEA_str()
    # Then we write all the oxidizers
    for ox in oxidizers:
      contents += ox.get_CEA_str()
    if self.inserts:
      contents += "insert "+" ".join(self.inserts)+"\n"
    if self.omits:
      contents += "omit "+" ".join(self.omits)+"\n"
      
    contents += "output {}trans\n".format("massf " if self.massf else "") # Output is mass fraction is "massf"
    contents += self.get_plt_string() # output plotting string, if any
    contents += "   plot cond\n"
    contents += "end\n"

    return contents
  
  def get_plt_string(self) -> str:
    return f"   plot {self.plt_keys:s}\n" # specify string formatting so None errors
  
  def _error_check_thermo_out_line(self, error_line:str):
    # Either raises an error or does nothing
    filePath = self.CEA.OUT_ERROR_FILE if self.CEA.dump_out_on_error else "(no output file dumped because of `dump_out_on_error` setting)"
    prepend = "CEA Failed to Run. "
    errors = { # Errors. Keys are error keys in CEA output. Values are a string to print which are prepended with a string and appended with the file name
    "FATAL": "FATAL error in input/output file: ",
    "LOW TEMPERATURE IMPLIES": "Reactants failed to react with given species \n  (likely need to add condensed phase product species to 'inserts')\n  Check output file for details: ",
    "CALCULATIONS STOPPED AFTER POINT": "No converged solution found for file: ",
    "CONVERGENCES FAILED": "CEA cannot determine correct condensed species. Try adding the correct ones with 'inserts'. Check output file: ",
    "DERIVATIVE MATRIX SINGULAR": "Encountered a singular matrix in output file: ",
    "ITERATIONS DID NOT SATISFY CONVERGENCE": "Failed to reach solution in reasonable number of iterations. Check output file: ",
    "ELECTRON BALANCE": "In ionic iteration, convergence criterion was not satisfied. Check output file: ",
    }
    
    for key, value in errors.items():
      if key in error_line:
        raise RuntimeError(prepend + value + filePath)
    
  def process_output(self, out_file_content:str, plt_file_content:str) -> OutputType:
    """
    Processes the output generated by CEA, taking the contents of the .out file and .plt file to generate output

    :param out_file_content: The contents of the .out file
    :param plt_file_content: The contents of the .plt file
    :raises RuntimeError: If the .plt file is empty
    :return: An Output with attributes defined by self.plt_keys
    """
    out = Output()
    if not plt_file_content: # If empty string
      raise RuntimeError("CEA Failed to Run. Plot file wasn't generated for " + self.filename)
    for line in plt_file_content.splitlines():
      new_line = line.split()
      if new_line[0] == '#':
        continue
      else:
        for key, value in zip(self._get_plt_keys().split(" "), new_line):
          out[key] = float(value)
    return out

  ## THESE MUST BE IMPLEMENTED BY SUBCLASSES ##
  def get_prefix_string(self) -> str:
    raise NotImplementedError()

@dataclass(kw_only=True)
class DetonationOutput(DictDataclass):
  prod_c: dict[str, float]
  massf: bool
  p: float
  t: float
  h: float
  rho: float
  son: float
  visc: float
  mw: float
  cp: float
  gammas: float
  gamma: float
  vel: float
  mach: float
  p_p1: float
  t_t1: float
  m_m1: float
  rho_rho1: float
  dLV_dLP_t: float
  dLV_dLT_p: float
  cond: float
  pran: float
  phi: float

class DetonationProblem(Problem[DetonationOutput]):
  problem_type = "det"
  plt_keys = "p t h mw cp gammas phi vel mach rho son cond pran"
  
  def get_prefix_string(self):
    to_ret = [
      str(self.problem_type),
      f"   p({self.pressure_units}) = {self.pressure:0.5f}",
      f"   {self.ratio_name} = {self.ratio_value:0.5f}",
    ]
    return "\n".join(to_ret) + "\n"

  def process_output(self, out_file_content:str, plt_file_content:str) -> DetonationOutput:
    out = Output()
    
    # We'll also open this file to get mass/mole fractions of all constituents and other properties
    # Ordered (in the file) list of terms that we're searching for
    search_terms = ["BURNED GAS", "DETONATION PARAMETERS", ("MASS FRACTIONS" if self.massf else "MOLE FRACTIONS")]
    search_i = 0
    in_search_area = False
    was_in_search_area = False # State variable so we increment search_i when it changes high to low
    passed_first_line = False # State variable so we ignore the empty line at the start of each section
    
    # The format of each section is "HEADING\n [empty line]\n [lines of data...]\n [empty line]
    #   So when we get a heading, we ignore the first line after, and keep going until the code for that section says to stop
    
    out.prod_c = Output()
    
    for line in out_file_content.splitlines():
      self._error_check_thermo_out_line(line)
      # If we are no longer in a search area, increment our search term to look for the next search area
      if was_in_search_area and not in_search_area:
        passed_first_line = False # reset this too
        search_i += 1
      was_in_search_area = in_search_area
      if search_i >= len(search_terms): # If we have run out of search terms, close the file
        break
    
      if search_terms[search_i] in line:
        in_search_area = True
      elif in_search_area and not passed_first_line:
        passed_first_line = True # Skip this iteration because its a blank line
      elif in_search_area and passed_first_line:
        # Lines should start with the first non-blank line
        # IMPORTANT: Each case should set in_search_area to False to end their section
        ### Combustion Chamber Products ###
        if search_i == 2:
          if not line.strip():
            in_search_area = False
            continue
          split = re.findall(r"\S+", line) # match by sections which are not whitespace
          key = split[0].lstrip("*") # remove any asterisks
          out.prod_c[key] = float(split[1])
        ### Detonation Products ###
        elif search_i==1:
          if not line.strip():
            in_search_area = False
            continue
          split = re.findall(r"\S+", line) # match by sections which are not whitespace
          key = split[0].replace("/","_").lower() # convert to usable format
          if key == "det": # This means we've gone past the p, t, m, rho / initial lines
            in_search_area = False
            continue # break out
          out[key] = float(split[1])
        ### Output gas products ###
        elif search_i == 0:
          if not line.strip(): # If this line is empty, skip it
            continue
          split = re.findall(r"\S+", line) # match by sections which are not whitespace
          key = split[0]
          if key == "(dLV/dLP)t":
            out.dLV_dLP_t = float(split[1])
          elif key == "(dLV/dLT)p":
            out.dLV_dLT_p = float(split[1])
          elif key == "VISC,MILLIPOISE": # Keep going until we find the viscosity line
            out.visc = float(split[1]) * 0.0001 # convert millipoise to Pa-s
            in_search_area = False
            continue
    
    out_plt = super().process_output(out_file_content, plt_file_content) # Process the plt file second so we can read errors listed in output file
    out.update(out_plt)
    
    # Also include the actual gamma and not just the isentropic gamma
    out.gamma = out.gammas*-out.dLV_dLP_t

    out.massf = self.massf
    
    return DetonationOutput(**out)

@dataclass(kw_only=True)
class HPOutput(DictDataclass):
  prod_c: dict[str, float]
  massf: bool
  p: float
  t: float
  h: float
  rho: float
  son: float
  visc: float
  mw: float
  cp: float
  gammas: float
  gamma: float
  dLV_dLP_t: float
  dLV_dLT_p: float
  cond: float
  pran: float
  phi: float

class HPProblem(Problem[HPOutput]):
  problem_type = "hp"
  plt_keys = "p t h mw cp gammas phi rho son cond pran"
  def get_prefix_string(self):
    to_ret = []
    to_ret.append("{}".format(self.problem_type))
    to_ret.append("   p({}) = {:0.5f}".format(self.pressure_units, self.pressure))
    to_ret.append("   {} = {:0.5f}".format(self.ratio_name, self.ratio_value))
    return "\n".join(to_ret) + "\n"
  
  def process_output(self, out_file_content:str, plt_file_content:str) -> HPOutput:
    out = Output()
    
    # We'll also open this file to get mass/mole fractions of all constituents and other properties
    # Ordered (in the file) list of terms that we're searching for
    search_terms = ["THERMODYNAMIC PROPERTIES", ("MASS FRACTIONS" if self.massf else "MOLE FRACTIONS")]
    search_i = 0
    in_search_area = False
    was_in_search_area = False # State variable so we increment search_i when it changes high to low
    passed_first_line = False # State variable so we ignore the empty line at the start of each section
    
    # The format of each section is "HEADING\n [empty line]\n [lines of data...]\n [empty line]
    #   So when we get a heading, we ignore the first line after, and keep going until the code for that section says to stop
    
    out.prod_c = Output()
    
    for line in out_file_content.splitlines():
      self._error_check_thermo_out_line(line)
      # If we are no longer in a search area, increment our search term to look for the next search area
      if was_in_search_area and not in_search_area:
        passed_first_line = False # reset this too
        search_i += 1
      was_in_search_area = in_search_area
      if search_i >= len(search_terms): # If we have run out of search terms, close the file
        break
    
      if search_terms[search_i] in line:
        in_search_area = True
      elif in_search_area and not passed_first_line:
        passed_first_line = True # Skip this iteration because its a blank line
      elif in_search_area and passed_first_line:
        # Lines should start with the first non-blank line
        # IMPORTANT: Each case should set in_search_area to False to end their section
        ### Combustion Chamber Products ###
        if search_i == 1:
          if not line.strip():
            in_search_area = False
            continue
          split = re.findall(r"\S+", line) # match by sections which are not whitespace
          key = split[0].lstrip("*") # remove any asterisks
          out.prod_c[key] = float(split[1])
        ### Output gas products ###
        elif search_i == 0:
          if not line.strip(): # If this line is empty, skip it
            continue
          split = re.findall(r"\S+", line) # match by sections which are not whitespace
          key = split[0]
          if key == "(dLV/dLP)t":
            out.dLV_dLP_t = float(split[1])
          elif key == "(dLV/dLT)p":
            out.dLV_dLT_p = float(split[1])
          elif key == "VISC,MILLIPOISE": # Keep going until we find the viscosity line
            out.visc = float(split[1]) * 0.0001 # convert millipoise to Pa-s
            in_search_area = False
            continue
    
    out_plt = super().process_output(out_file_content, plt_file_content) # Process the plt file second so we can read errors listed in output file
    out.update(out_plt)
    
    # Also include the actual gamma and not just the isentropic gamma
    out.gamma = out.gammas*-out.dLV_dLP_t

    out.massf = self.massf
    
    return HPOutput(**out)

class TPMaterial(Material):
  def __init__(self, parent: Material):
    try:
      p = parent
      self.parent = p
      self.output_type = self.parent.output_type
      self.is_fuel = self.parent.is_fuel
      super().__init__(p._name, p._temp, p._wt, p._mols, p._chemical_representation)
    except AttributeError:
      raise TypeError("TPMaterial must be derived from Material, but a '{}' object was supplied".format(type(parent)))

  def get_CEA_name_wt(self, wt_or_mol: str, amount: float) -> str:
    """ Like the superclass, but we don't supply a temperature per-input for TP problems """
    return "   {}={} {}={:0.5f}".format(self.output_type, self.name, wt_or_mol, amount)

class TPProblem(HPProblem):
  problem_type = "tp"
  
  def __init__(self, *args, 
    temperature=298, # Reaction temperature, in Kelvin
    temperature_units="k", # units for temperature
    **kwargs
  ):
    super().__init__(*args, **kwargs)
      
    self.temperature = temperature
    self._temperature_units = None
    self.set_temperature_units(temperature_units)
  
  def set_temperature(self, temperature: float): self.temperature = temperature
  
  def set_temperature_units(self, temperature_units: str):
    temperature_units = temperature_units.lower()  # We only have a few options for pressure units
    if temperature_units not in OPTIONS_TEMP_UNITS:
      raise ValueError("pressure unit must be in " + ", ".join(OPTIONS_TEMP_UNITS))
    self._temperature_units = temperature_units
  
  @property
  def temperature_units(self): return self._temperature_units
  @temperature_units.setter
  def _temperature_units_setter(self, value: str): return self.set_temperature_units(value)


  def get_prefix_string(self):
    to_ret = []
    to_ret.append("{}".format(self.problem_type))
    to_ret.append("   p({}) = {:0.5f}".format(self.pressure_units, self.pressure))
    to_ret.append("   t({}) = {:0.5f}".format(self._temperature_units, self.temperature))
    to_ret.append("   {} = {:0.5f}".format(self.ratio_name, self.ratio_value))
    return "\n".join(to_ret) + "\n"
  
  def make_input_file(self, material_list: list[Material]): # chamber conditions and materials list
    # Here we convert all the materials to TPMaterials before calling the main function, so that each one doesn't report temperature
    new_material_list = [TPMaterial(material) for material in material_list]
    return super().make_input_file(new_material_list)

def make_area_property(name: str) -> property:
  """ Helper for RocketProblem objects to make properties for all the different area ratios """
  def getter(self: "RocketProblem"): 
    f"""
    Get/set the "{name}" nozzle ratio property for this RocketProblem
    If the current nozzle ratio property is not "{name}" and you try to get the value, a ValueError will be raised
    """
    if self.nozzle_ratio_name == name:
      return self.ratio_value
    raise ValueError(f'Currently set area ratio name was "{self.nozzle_ratio_name}", not {name}')
  
  def setter(self: "RocketProblem", value: float):
    return self._set_area_ratio(name, value)
  
  return property(getter, setter)

@dataclass(kw_only=True)
class BaseRocketOutput(DictDataclass):
  """
  Data representation for rocket engine properties.

  :param prod_c: Chamber products
  :param prod_t: Throat products
  :param prod_e: Exit products
  :param p: Pressure, bar
  :param t_p: Throat pressure
  :param c_p: Chamber pressure
  :param t: Temperature, Kelvin
  :param t_t: Throat temperature
  :param c_t: Chamber temperature
  :param h: Enthalpy, kJ/kg
  :param t_h: Throat enthalpy
  :param c_h: Chamber enthalpy
  :param rho: Density, kg/m^3
  :param t_rho: Throat density
  :param c_rho: Chamber density
  :param son: Sonic velocity, m/s
  :param t_son: Throat sonic velocity
  :param c_son: Chamber sonic velocity
  :param visc: Burned gas viscosity, Pascal-Seconds
  :param t_visc: Throat viscosity
  :param c_visc: Chamber viscosity
  :param cond: Burned gas thermal conductivity, W/(m K)
  :param t_cond: Throat thermal conductivity
  :param c_cond: Chamber thermal conductivity
  :param pran: Burned gas Prandtl number
  :param t_pran: Throat Prandtl number
  :param c_pran: Chamber Prandtl number
  :param mw: Molecular weight of all products, kg/kmol
  :param t_mw: Throat molecular weight
  :param c_mw: Chamber molecular weight
  :param m: molecular weight calculated as the weight of all products divided by the number of gaseous moles (same as mw if no condensed phases as mw=m), kg/kmol
  :param t_m: Throat molecular weight
  :param c_m: Chamber molecular weight
  :param condensed: True if condensed phase products exist
  :param t_condensed: Throat condensed phase existence
  :param c_condensed: Chamber condensed phase existence
  :param cp: Constant-pressure specific heat capacity, kJ/(kg*K)
  :param t_cp: Throat cp
  :param c_cp: Chamber cp
  :param gammas: Isentropic ratio of specific heats
  :param t_gammas: Throat gammas
  :param c_gammas: Chamber gammas
  :param gamma: Real ratio of specific heats
  :param t_gamma: Throat gamma
  :param c_gamma: Chamber gamma
  :param isp: Ideal ISP (ambient pressure = exit pressure), s
  :param t_isp: Throat ISP
  :param ivac: Vacuum ISP, s
  :param t_ivac: Throat vacuum ISP
  :param cf: Ideally expanded thrust coefficient
  :param t_cf: Throat CF
  :param dLV_dLP_t: (dLV/dLP)t
  :param t_dLV_dLP_t: Throat dLV/dLP
  :param c_dLV_dLP_t: Chamber dLV/dLP
  :param dLV_dLT_p: (dLV/dLT)p
  :param t_dLV_dLT_p: Throat dLV/dLT
  :param c_dLV_dLT_p: Chamber dLV/dLT
  :param cstar: Characteristic velocity in chamber, m/s
  :param mach: Mach number at exhaust
  :param o_f: Oxidizer/Fuel weight ratio
  :param phi: Weight-based equivalence ratio of oxidizer/fuel
  :param ae: Ratio of area at exit to area at throat
  :param t_ae: Throat area ratio (always 1)
  :param pip: Ratio of pressure in chamber to pressure at exit
  :param t_pip: Ratio of pressure in chamber to pressure at throat
  """
  prod_c: Output
  massf: bool
  p: float
  t_p: float
  c_p: float
  t: float
  t_t: float
  c_t: float
  h: float
  t_h: float
  c_h: float
  rho: float
  t_rho: float
  c_rho: float
  son: float
  t_son: float
  c_son: float
  visc: float
  t_visc: float
  c_visc: float
  cond: float
  t_cond: float
  c_cond: float
  pran: float
  t_pran: float
  c_pran: float
  mw: float
  t_mw: float
  c_mw: float
  m: float
  t_m: float
  c_m: float
  condensed: bool
  t_condensed: bool
  c_condensed: bool
  cp: float
  t_cp: float
  c_cp: float
  gammas: float
  t_gammas: float
  c_gammas: float
  gamma: float
  t_gamma: float
  c_gamma: float
  isp: float
  t_isp: float
  ivac: float
  t_ivac: float
  cf: float
  t_cf: float
  
  cstar: float
  mach: float
  o_f: float
  phi: float
  ae: float
  t_ae: float
  pip: float
  t_pip: float

@dataclass(kw_only=True)
class FrozenRocketOutput(BaseRocketOutput):
  """
  Frozen flows have a few extra elements because the equilibrium values are *also* calculated and they may be important to someone

  :param cond_eq: Burned gas thermal conductivity, W/(m K) in equilibrium flow at this station
  :param t_cond_eq: Throat equilibrium conductivity
  :param c_cond_eq: Chamber equilibrium conductivity
  """
  cond_eq: float
  t_cond_eq: float
  c_cond_eq: float

  '''
  As of the time of writing, there is a limit of 19 plt output keys, and so we can't add equilibrium prandtl number to the plt.
  I am keeping these docstrings for posterity
  
  :param pran_eq: Burned gas Prandtl number in equilibrium flow at this station
  :param t_pran_eq: Throat equilibrium Prandtl number
  :param c_pran_eq: Chamber equilibrium Prandtl number
  '''
  # pran_eq: float
  # t_pran_eq: float
  # c_pran_eq: float

@dataclass(kw_only=True)
class RocketOutput(BaseRocketOutput):
  """
  Equilibrium problems have access to several other attributes, not present in Frozen flow

  :param prod_t: Throat products
  :param prod_e: Exit products
  :param dLV_dLP_t: (dLV/dLP)t
  :param t_dLV_dLP_t: Throat dLV/dLP
  :param c_dLV_dLP_t: Chamber dLV/dLP
  :param dLV_dLT_p: (dLV/dLT)p
  :param t_dLV_dLT_p: Throat dLV/dLT
  :param c_dLV_dLT_p: Chamber dLV/dLT
  """
  prod_t: Output
  prod_e: Output
  dLV_dLP_t: float
  t_dLV_dLP_t: float
  c_dLV_dLP_t: float
  dLV_dLT_p: float
  t_dLV_dLT_p: float
  c_dLV_dLT_p: float

@dataclass(kw_only=True)
class FiniteAreaCombustorRocketOutput(RocketOutput):
  """
  FiniteAreaCombustor problems have access to products at the combustor end

  :param prod_f: Products at the end of the combustor
  :param f_dLV_dLP_f: (dLV/dLP)f at combustor end
  :param f_dLV_dLT_p: (dLV/dLT)p at combustor end
  :param f_visc: Burned gas viscosity at combustor end, Pascal-Seconds
  :param f_isp: Ideal ISP at combustor end (ambient pressure = exit pressure), s
  :param f_ivac: Vacuum ISP at combustor end, s
  :param f_cf: Ideally expanded thrust coefficient at combustor end
  :param f_p: Pressure at combustor end, bar
  :param f_t: Temperature at combustor end, Kelvin
  :param f_m: Molecular weight considering condensed phases at combustor end, kg/kmol
  :param f_mw: Molecular weight of all products at combustor end, kg/kmol
  :param f_condensed: True if condensed phase products exist at combustor end
  :param f_cp: Constant-pressure specific heat capacity at combustor end, kJ/(kg*K)
  :param f_gammas: Isentropic ratio of specific heats at combustor end
  :param f_rho: Density at combustor end, kg/m^3
  :param f_son: Sonic velocity at combustor end, m/s
  :param f_h: Enthalpy at combustor end, kJ/kg
  :param f_cond: Burned gas thermal conductivity at combustor end, W/(m K)
  :param f_pran: Burned gas Prandtl number at combustor end
  :param f_ae: Ratio of area at combustor end to area at throat
  :param f_pip: Ratio of pressure in chamber to pressure at combustor end
  """
  prod_f: Output
  f_dLV_dLP_f: float
  f_dLV_dLT_p: float
  f_visc: float
  f_isp: float
  f_ivac: float
  f_cf: float
  f_p: float
  f_t: float
  f_m: float
  f_mw: float
  f_condensed: bool
  f_cp: float
  f_gammas: float
  f_rho: float
  f_son: float
  f_h: float
  f_cond: float
  f_pran: float
  f_ae: float
  f_pip: float

class RocketProblem(Problem[RocketOutput|FrozenRocketOutput|FiniteAreaCombustorRocketOutput]):
  problem_type = "rocket"
  plt_keys = "p t isp ivac m mw cp gam o/f cf rho son mach phi h cond pran ae pip"
    
  def __init__(self, *args, sup: float|None=None, sub: float|None=None, ae_at: float|None=None, pip: float|None=None, 
               analysis_type:str="equilibrium", fac_ac:float|None=None, fac_ma:float|None=None, nfz: int|None=None,
               custom_nfz: float|None=None,
               **kwargs
               ):
    super().__init__(*args, **kwargs)
    
    if ae_at:
      sup = ae_at
    if not sup and not sub and not pip: # Default if nothing is specified
      sup = 1
    if sup and sub:
      raise ValueError("Can only specify supersonic or subsonic area ratio, not both")
    if pip and (sup or sub):
      raise ValueError("Can only specify area ratio or pressure ratio, not both")
    if fac_ac and fac_ma:
      raise ValueError("Can only specify fac ac/at or fac ma, not both")

     # ensure case because we check for frozen by literal
    self.nozzle_ratio_name = "pip" if pip else ("sup" if sup else "sub") # Can only specify one ratio, pressure or supersonic or subsonic area ratio
    self.nozzle_ratio_value = [pip if pip else self._validate_ratio(sup if sup else sub), ] # 1-element array
    self.fac_type = None
    self.fac_value = None
    self.analysis_type: str  # Specify before calling
    self.set_analysis_type(analysis_type, nfz, custom_nfz=custom_nfz)
    if fac_ac: self.set_fac_ac(fac_ac)
    if fac_ma: self.set_fac_ma(fac_ma)
    
  
  
  @staticmethod
  def _validate_ratio(ratio):
    if ratio < 1:
      raise ValueError("Supersonic/Subsonic area ratio must always be >= 1")
    return ratio
  def _set_area_ratio(self, ratio_name: str, value: float):
    self.nozzle_ratio_name = ratio_name
    if ratio_name in {"sup", "sub"}: value = self._validate_ratio(value)
    elif ratio_name == "pip": pass
    else: raise ValueError(f"Unknown ratio name '{ratio_name}'")
    self.nozzle_ratio_value[-1] = value

  def set_sup(self, sup: float): self._set_area_ratio("sup", sup)
  def set_ae_at(self, ae_at): self.set_sup(ae_at)
  def set_sub(self, sub: float): self._set_area_ratio("sub", sub)
  def set_pip(self, pip: float): self._set_area_ratio("pip", pip)

  # Make properties for each type of area ratio so you can do `problem.sup = 5` or similar
  sup = make_area_property("sup")
  ae_at = make_area_property("sup") # Also sup
  sub = make_area_property("sub")
  pip = make_area_property("pip")
  
  def set_analysis_type(self, analysis_type: str, nfz=None, custom_nfz=None):
    analysis_type = analysis_type.lower()
    if nfz is not None and custom_nfz is not None:
      raise ValueError("cannot specify both nfz and custom_nfz. nfz= will be set automatically if custom_nfz is set")
    if "equilibrium" in analysis_type and "frozen" in analysis_type:
      raise ValueError("Rocket_Problem does not support combined equilibrium-frozen calculations")
    if self.fac_type and "frozen" in analysis_type:
      raise ValueError("Rocket_Problem does not support combined finite area combustor and frozen calculations")

    # Anytime we change analysis type, reset the nozzle ratio to be only the last value
    self.nozzle_ratio_value = self.nozzle_ratio_value[-1:]  # List with only the last value
    if "frozen" in analysis_type:
      if custom_nfz is not None:
        nfz = 3 # 3 is the column for the custom nfz point
        if custom_nfz > self.nozzle_ratio_value[-1]:
          raise ValueError("nozzle ratio for the frozen point must be lower than the outlet nozzle ratio")
        self.nozzle_ratio_value.insert(0, custom_nfz)  # Add another column at the point we want to freeze

      if isinstance(nfz, int):
        if "nfz" in analysis_type:  # If they already had "nfz" in analysis_type, let them know
          raise ValueError("already had 'nfz=' in analysis type. Cannot specify multiple frozen locations")
        analysis_type = analysis_type + " nfz={}".format(nfz)
      # Modify plot keys to use frozen values
      to_put = type(self)._get_plt_keys() # Get copy from class instance
      for orig in ("isp", "ivac", "cf", "mach", "cond", "pran"):
        to_put = to_put.replace(" "+orig, " "+orig+"fz", 1) # Add 'fz' to each of these properties
      to_put += " cond" # Add in the original equilibrium values for conductivity as well, it can be important. I would like to also add "pran" but there is a limit on columns in the .plt file
      self.plt_keys = to_put
    else:
      self.plt_keys = type(self).plt_keys # Get new copy from class instance
    self.analysis_type = analysis_type
  
  def _set_fac_check(self): # Simple check to make sure we don't set this.
    if self.analysis_type is not None and "frozen" in self.analysis_type:
      raise ValueError("Rocket_Problem does not support combined finite area combustor and frozen calculations")
  def set_fac_ac(self, fac): self._set_fac_check(); self.fac_type = "ac/at"; self.fac_value = fac
  def set_fac_ma(self, fac): self._set_fac_check(); self.fac_type = "mdot"; self.fac_value = fac
  def unset_fac(self): self.fac_type = None; self.fac_value = None
  
  def get_prefix_string(self):
    to_ret = list()
    to_ret.append("{} {}".format(self.problem_type, self.analysis_type))
    to_ret.append("   p({}) = {:0.5f}".format(self.pressure_units, self.pressure))
    to_ret.append("   {} = {:0.5f}".format(self.ratio_name, self.ratio_value))
    # enable using multiple area ratio for custom frozen point
    to_ret.append("   {} = ".format(self.nozzle_ratio_name) + ",".join("{:0.5f}".format(rat) for rat in self.nozzle_ratio_value))
    if self.fac_type:
      to_ret.append("   fac {} = {:0.5f}".format(self.fac_type, self.fac_value))
    return "\n".join(to_ret) + "\n"
  
  def process_output(self, out_file_content: str, plt_file_content: str) -> RocketOutput|FrozenRocketOutput|FiniteAreaCombustorRocketOutput:
    out = Output()

    is_frozen = "frozen" in self.analysis_type
    
    # We'll also open this file to get mass/mole fractions of all constituents and other properties
    # Ordered (in the file) list of terms that we're searching for
    search_terms = ["CHAMBER", "PERFORMANCE PARAMETERS", ("MASS FRACTIONS" if self.massf else "MOLE FRACTIONS")]
    search_i = 0
    in_search_area = False
    was_in_search_area = False # State variable so we increment search_i when it changes high to low
    passed_first_line = False # State variable so we ignore the empty line at the start of each section
    
    # The format of each section is "HEADING\n [empty line]\n [lines of data...]\n [empty line]
    #   So when we get a heading, we ignore the first line after, and keep going until the code for that section says to stop
    
    has_fac = bool(self.fac_type) # If has finite area combustor
    ch_th = 3 if has_fac else 2
    end_col = ch_th + len(self.nozzle_ratio_value)
    
    out.prod_c = Output()
    out.prod_t = Output()
    out.prod_e = Output()
    if has_fac:
      out.prod_f = Output()
    
    # Float map, mapping columns because the mapping isn't always 1,2,3
    def flMap(split_line, key):
      if has_fac: # If finite area combustor, we have more columns
        mapping = {"c": 1, "f": 2, "t": 3, "e": end_col}
      else: # Otherwise, normal
        mapping = {"c": 1, "t": 2, "e": end_col}
      return float(split_line[mapping[key]])
    
    for line in out_file_content.splitlines():
        self._error_check_thermo_out_line(line)
        # If we are no longer in a search area, increment our search term to look for the next search area
        if was_in_search_area and not in_search_area:
          passed_first_line = False # reset this too
          search_i += 1
        was_in_search_area = in_search_area
        if search_i >= len(search_terms): # If we have run out of search terms, close the file
          break
      
        if search_terms[search_i] in line:
          in_search_area = True
        elif in_search_area and not passed_first_line: # After every section title is a blank line
          passed_first_line = True # Skip this iteration because its a blank line
        elif in_search_area and passed_first_line:
          # Lines should start with the first non-blank line
          # IMPORTANT: Each case should set in_search_area to False to end their section
          ### Combustion Chamber Products ###
          if search_i == 2:
            if not line.strip():
              in_search_area = False
              continue
            split = re.findall(r"\S+", line) # match by sections which are not whitespace
            # So apparently in frozen calculations, the products come in multiple columns and only chamber
            if "frozen" in self.analysis_type:
              # Split them into pairs 2-by-2
              kv_pairs = [(split[i], split[i+1]) for i in range(0, len(split), 2)]
              for key, value in kv_pairs:
                key = key.lstrip("*") # remove any asterisks
                out.prod_c[key] = float(value)
            else:
              key = split[0].lstrip("*") # remove any asterisks
              out.prod_c[key] = flMap(split, 'c')
              out.prod_t[key] = flMap(split, 't')
              out.prod_e[key] = flMap(split, 'e')
              if has_fac: out.prod_f[key] = flMap(split, 'f')
              
          ### Output gas products ###
          elif search_i == 0:
            if not line.strip(): # If this line is empty, skip it
              continue
            split = re.findall(r"\S+", line) # match by sections which are not whitespace
            key = split[0]
            if key == "(dLV/dLP)t": # Will not exist in frozen
              out.c_dLV_dLP_t = flMap(split, 'c')
              out.t_dLV_dLP_t = flMap(split, 't')
              out.dLV_dLP_t = flMap(split, 'e')
              if has_fac: out.f_dLV_dLP_f = flMap(split, 'f')
            elif key == "(dLV/dLT)p":# Will not exist in frozen
              out.c_dLV_dLT_p = flMap(split, 'c')
              out.t_dLV_dLT_p = flMap(split, 't')
              out.dLV_dLT_p = flMap(split, 'e')
              if has_fac: out.f_dLV_dLT_p = flMap(split, 'f')
            elif key == "VISC,MILLIPOISE": # Keep going until we find the viscosity line
              out.c_visc = flMap(split, 'c') * 0.0001 # convert millipoise to Pa-s
              out.t_visc = flMap(split, 't') * 0.0001 # convert millipoise to Pa-s
              out.visc = flMap(split, 'e') * 0.0001 # convert millipoise to Pa-s
              if has_fac: out.f_visc = flMap(split, 'f')
              in_search_area = False
              continue
          elif search_i == 1:
            if not line.strip(): # If this line is empty, skip it
              in_search_area = False
              continue
            split = re.findall(r"\S+", line) # match by sections which are not whitespace
            key = split[0]
            if key == "CSTAR,":
              out.cstar = float(split[2])
            if "frozen" in self.analysis_type: # These don't show up in the plot file...
              if key == "CF":
                out.t_cf = float(split[1])
                out.cf = float(split[2])
              elif key == "Ivac,":
                out.t_ivac = float(split[2])/9.81
                out.ivac = float(split[3])/9.81
              elif key == "Isp,":
                out.t_isp = float(split[2])/9.81
                out.isp = float(split[3])/9.81
    
    if not plt_file_content:
      raise RuntimeError("CEA Failed to Run. Plot file wasn't generated")
    
    chamber_count = 1
    if has_fac:
      fac_count = 2
      throat_count = 3
      exit_count = 4
    else:
      fac_count = -1
      throat_count = 2
      exit_count = end_col
    
    for counter, line in enumerate(plt_file_content.splitlines()):
          new_line = line.split()
          if counter == 0:
            continue
          elif  counter == chamber_count: # chamber properties
            # Values common to all
            out.o_f = float(new_line[8])
            out.phi = float(new_line[13])
            
            # Section Properties 
            out.c_p = float(new_line[0])
            out.c_t = float(new_line[1])
            out.c_m = float(new_line[4])
            out.c_mw = float(new_line[5])
            if out.c_m == 0: # If no condensed phase products, this becomes 0 because its the same as mw
              out.c_condensed = False
              out.c_m = out.c_mw
            else:
              out.c_condensed = True # has condensed phase products
            out.c_cp = float(new_line[6])
            out.c_gammas = float(new_line[7])
            out.c_rho = float(new_line[10])
            out.c_son = float(new_line[11])
            out.c_h = float(new_line[14])
            out.c_cond = float(new_line[15])/10 # Convert mW/cm K to W/m K
            out.c_pran = float(new_line[16])
            if is_frozen: # Added on to end
              out.c_cond_eq = float(new_line[19])
          elif counter == fac_count: # finite area combustor properties
            # Throat/Nozzle Values
            out.f_isp = float(new_line[2])/9.81
            out.f_ivac = float(new_line[3])/9.81
            out.f_cf = float(new_line[9])

            # Section Properties  
            out.f_p = float(new_line[0])
            out.f_t = float(new_line[1])
            out.f_m = float(new_line[4])
            out.f_mw = float(new_line[5])
            if out.f_m == 0: # If no condensed phase products, this becomes 0 because its the same as mw
              out.f_condensed = False
              out.f_m = out.f_mw
            else:
              out.f_condensed = True # has condensed phase products
            out.f_cp = float(new_line[6])
            out.f_gammas = float(new_line[7])
            out.f_rho = float(new_line[10])
            out.f_son = float(new_line[11])
            out.f_h = float(new_line[14])
            out.f_cond = float(new_line[15])/10
            out.f_pran = float(new_line[16])
            out.f_ae = float(new_line[17])
            out.f_pip = float(new_line[18])
          elif counter == throat_count: # throat properties
            # Throat/Nozzle Values
            out.t_isp = float(new_line[2])/9.81
            out.t_ivac = float(new_line[3])/9.81
            out.t_cf = float(new_line[9])
            
            # Section Properties  
            out.t_p = float(new_line[0])
            out.t_t = float(new_line[1])
            out.t_m = float(new_line[4])
            out.t_mw = float(new_line[5])
            if out.t_m == 0: # If no condensed phase products, this becomes 0 because its the same as mw
              out.t_condensed = False
              out.t_m = out.t_mw
            else:
              out.t_condensed = True # has condensed phase products
            out.t_cp = float(new_line[6])
            out.t_gammas = float(new_line[7])
            out.t_rho = float(new_line[10])
            out.t_son = float(new_line[11])
            out.t_h = float(new_line[14])
            out.t_cond = float(new_line[15])/10
            out.t_pran = float(new_line[16])
            out.t_ae = float(new_line[17])
            out.t_pip = float(new_line[18])
            if is_frozen: # Added on to end
              out.t_cond_eq = float(new_line[19])
          elif counter == exit_count: # nozzle exit properties
            # Exit-only values
            out.mach = float(new_line[12])
            
            # Throat/Nozzle Values
            out.isp = float(new_line[2])/9.81
            out.ivac = float(new_line[3])/9.81
            out.cf = float(new_line[9])

            # Section Properties  
            out.p = float(new_line[0])
            out.t = float(new_line[1])
            out.m = float(new_line[4])
            out.mw = float(new_line[5])
            if out.m == 0: # If no condensed phase products, this becomes 0 because its the same as mw
              out.condensed = False
              out.m = out.mw
            else:
              out.condensed = True # has condensed phase products
            out.cp = float(new_line[6])
            out.gammas = float(new_line[7])
            out.rho = float(new_line[10])
            out.son = float(new_line[11])
            out.h = float(new_line[14])
            out.cond = float(new_line[15])/10
            out.pran = float(new_line[16])
            out.ae = float(new_line[17])
            out.pip = float(new_line[18])
            if is_frozen: # Added on to end
              out.cond_eq = float(new_line[19])

    if is_frozen:
      # We don't get mole fractions in the other positions for frozen
      del out["prod_t"]
      del out["prod_e"]
      # In frozen the normal equals isentropic
      out.gamma, out.c_gamma, out.t_gamma = out.gammas, out.c_gammas, out.t_gammas
    else:
      # Also include the actual gamma and not just the isentropic gamma
      out.gamma = out.gammas*-out.dLV_dLP_t
      out.c_gamma = out.c_gammas*-out.c_dLV_dLP_t
      out.t_gamma = out.t_gammas*-out.t_dLV_dLP_t

    out.massf = self.massf

    if is_frozen:
      return FrozenRocketOutput(**out)
    elif has_fac:
      return FiniteAreaCombustorRocketOutput(**out)
    else:
      return RocketOutput(**out)
