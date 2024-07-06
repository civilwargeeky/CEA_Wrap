import subprocess, re
import logging as _logging
import platform
import warnings
from typing import List, Union
from .utils import _get_data_file, cleanup_package_install, move_file_if_changed, Output
from .thermo_lib import ThermoInterface
logging = _logging.getLogger(__name__)

_BASE_CEA = "FCEA2.exe" if platform.system() == "Windows" else "FCEA2"

# The first time we install from source, we need to move our files from the ".assets" directory to our data directory
# This function checks if the assets directory still exists, and if so, will try to copy files to it.
cleanup_package_install()

CEA_LOCATION = _get_data_file(_BASE_CEA) # Get data file after cleanup because location may be different depending on availability
OPTIONS_TEMP_UNITS = ["k", "r", "c", "f"] # options for temperature units
OPTIONS_PRES_UNITS = ["bar", "atm", "psi", "mmh"] # options for pressure units
OPTIONS_DENS_UNITS = ["kg","g"] # options for Density units


def reload_thermo_lib():
  try:
    for file in ["thermo.lib", "trans.lib"]:
      pack_file = _get_data_file(file)
      move_file_if_changed(file, pack_file)
  except PermissionError as e:
    logging.error("---- Error! Attempted to copy thermo.lib and trans.lib into current directory but failed ----")
    logging.error("---- Is your current directory system32 or another protected directory? ----")
    raise e from None
  # Load our interface to all the ThermoMaterials
  ThermoInterface.load()
reload_thermo_lib() # Call this immediately

class Material:
  output_type = None # MUST BE SPECIFIED IN SUBCLASSES. The string we put before the material when writing .inp files
  is_fuel = None # MUST BE SPECIFIED IN SUBCLASSES. Assume all materials are either fuels or oxidizers. Must be overriden
  # If true, on initialization, we check if the material exists in the supplied thermo_spg input file
  #   also, will add .ref member to all objects for the thermo input reference
  check_against_thermo_inp = True 
  
  def __init__(self, name, temp:float=298.15, wt_percent:float=None, mols:float=None, chemical_composition:str=None, hf:float=None, hf_type="mol"):
    if wt_percent is None and mols is None: # If neither is specified, user probably doesn't care
      wt_percent = 100
    
    self.name = name # Name to input into CEA
    self.temp = temp # temperature in K, defaults to 298K
    self.wt_percent = wt_percent # Should be set when actually running through
    self.mols = mols # Mole ratio for this material
    self.chemical_composition = chemical_composition # chemical composition such as "LI 1 B 1 H 4" for LiBH4. If defined, will not use CEA default values
    self.hf = hf # Enthalpy of formation, if needs to be defined.
    if (chemical_composition == None) ^ (hf == None): # XOR is true if exactly one is true
      raise ValueError("Elements entered with exploded chemical formula or hf must have both")
    
    if self.check_against_thermo_inp and not chemical_composition: # If they don't override thermo.lib data and we check for missing elements
      if name not in ThermoInterface:
        close_matches = ThermoInterface.get_close_matches(name)
        raise ValueError(f"specified element '{name}' does not exist in package thermo library\n" +
                         f"Change name or set {__package__}.Material.check_against_thermo_inp to False\n"+
                         f"Material '{name}' {len(close_matches)} closest matches: \"" + '", "'.join(close_matches)+'"')
      
      self.ref = ThermoInterface[name] # Get the ThermoMaterial and store it in ref
    
      if not ThermoInterface[name].defined_at(temp):
        ranges = ThermoInterface[name].temp_ranges
        string = f"specified material '{name}' does not exist at temperature {temp:0.2f} " +\
                  " (min, max) = ({:.2f}, {:0.2f})".format(min(ranges)[0], max(ranges)[1]) # min and max sort by 1st element of tuples, assume tuple with max 1st element is max total.
        raise ValueError(string)
    
    if wt_percent and mols:
      raise TypeError("Material cannot have both wt_percent and mols specified")
  
  def set_wt_percent(self, wt_percent:float):
    self.mols = None # Can only have one
    self.wt_percent = wt_percent
  
  def set_mols(self, mols:float):
    self.wt_percent = None # Can only have one
    self.mols = mols
    
  def set_temp(self, temp:float):
    if self.check_against_thermo_inp and not ThermoInterface[self.name].defined_at(temp):
      raise ValueError(f"specified material '{self.name}' does not exist at temperature {temp:0.2f}")
    self.temp=temp
    
  def is_mols(self) -> bool: # Helper function for a material being in wt_percent or mols
    return self.mols is not None
    
  def is_wt_percent(self) -> bool: # Helper function for a material being in wt_percent or mols
    return self.wt_percent is not None

  def get_CEA_name_wt(self, wt_or_mol: str, amount: float) -> str:
    """ Gives the actual line put into .inp files, minus potential hf strings or chemical composition """
    return "   {}={} {}={:0.5f}  t,k={:0.5f}".format(self.output_type, self.name, wt_or_mol, amount, self.temp)

  def get_CEA_str(self) -> str:
    # Specify whether to use the str/val for weight or mols.
    name, ratio = "wt" if self.wt_percent is not None else "mol", self.wt_percent if self.wt_percent is not None else self.mols
    if ratio < 0:
      raise ValueError("Cannot have {} with <0 {} percent!".format(self.name, name))
    elif ratio == 0: # This allows us to not include materials that are set to 0 percent
      return ""
    else:
      string = self.get_CEA_name_wt(name, ratio)
      if self.hf:
        string += " h,kj/mol={:0.5f}".format(self.hf)
      if self.chemical_composition:
        string += " " + self.chemical_composition
      return string + " \n"
      
  def __str__(self):
    return self.name
      
class Fuel(Material):
  output_type = "fuel"
  is_fuel = True

F = Fuel # Alias

class Oxidizer(Material):
  output_type = "ox"
  is_fuel = False

O = Oxidizer # Alias

## Plan: Can also have composite Fuel/Oxidizer made up of percentages of other components

def run_cea_backend(filename:str):
  ret = subprocess.run(CEA_LOCATION, input=filename+"\n", text=True, stdout=subprocess.DEVNULL)
  if ret.returncode != 0:
    logging.error(ret)
    raise RuntimeError("Running CEA failed with errors")
  return ret
  
# My idea is that we could have different types of problems with similar methods for like get_prefix_string and things
class Problem:
  ### NOTE : ALL PROBLEM SUBCLASSES MUST SPECIFY PROBLEM TYPE AND PLOT KEYS ###
  problem_type = None
  plt_keys = None # plt_keys should be a space-separated string of items to go into the "plt" command of FCEA2
  
  # If true, on initialization, we check if all inserts and omits exist in the supplied thermo_spg input file
  check_against_thermo_inp = True 
  
  _ratio_options = ["p_f", "f_o", "o_f", "phi", "r_eq"] # Possible function arguments
  _ratio_CEA =     ["%f",  "f/o", "o/f", "phi", "r"] # Values to put into CEA
  def _set_fuel_ratio(self, **kwargs):
    # This function will go through kwargs and find which ratio we should use
    found_one = False
    for option, CEA in zip(self._ratio_options, self._ratio_CEA):
      if option in kwargs:
        if found_one: # If we already found one, we can't use another
          raise TypeError("Can only specify one ratio at once (%f, o/f, phi, etc.)")
        self.ratio_name = CEA # put in the string to place into CEA
        self.ratio_value = kwargs[option] # Then get the float
        found_one = True
  
  def _format_input_list(self, inputList):
    # Transmutes the inputList to a space-separated list of names
    # checks each element to ensure it is a valid material
    if inputList is None:
      return None
    if isinstance(inputList, str):
      inputList = inputList.split()
    inputList = list(map(str, inputList)) # make any Materials supplied into strings
    if self.check_against_thermo_inp:
      for material in inputList:
        if material not in ThermoInterface:
          close_matches = ThermoInterface.get_close_matches(material)
          raise ValueError(f"specified element '{material}' does not exist in package thermo library\n" +
                           f"Change name or set {__package__}.{__class__}.check_against_thermo_inp to False\n"+
                           f"Material '{material}' {len(close_matches)} closest matches: \"" + '", "'.join(close_matches)+'"')
    return inputList
    
  # All arguments must be specified by keyword
  def __init__(self, *,
               pressure: float=1000,  # Chamber/operation pressure
               materials: List[Material]=None,  # List of Material objects
               massf: bool=False,  # mass fractions or mol fractions in output
               filename: str="my_output",  # The file to be used for .inp/.out/.plt files
               pressure_units: str="psi",  # units for pressure
               density_units: str="kg", # units for density
               inserts: Union[str, List[Union[str, Material]]]=None,  # space-separated string or list of inserts
               omits:   Union[str, List[Union[str, Material]]]=None,  # space-separated string or list of omits
               **kwargs
               ):
      
    self.massf = massf
    self.materials = materials
    self.pressure = pressure
    self.pressure_units = None
    self.set_pressure_units(pressure_units)
    self.filename = None
    self.set_filename(filename)
    
    # Check against thermo file for these to prevent errors
    self.inserts = self._format_input_list(inserts)
    self.omits   = self._format_input_list(omits)
    
    self.ratio_name = None
    self.ratio_value = None
    self._set_fuel_ratio(**kwargs)
    
    diff = set(kwargs).difference(set(self._ratio_options))
    if diff: # If there are any kwargs keys that aren't in _ratio_options keys, we should error
      raise TypeError(self.__class__.__name__+"() got an unexpected keyword argument(s): " + ",".join(diff))
    
  def set_p_f(self, p_f: float): self._set_fuel_ratio(p_f=p_f)
  def set_f_o(self, f_o: float): self._set_fuel_ratio(f_o=f_o)
  def set_o_f(self, o_f: float): self._set_fuel_ratio(o_f=o_f)
  def set_phi(self, phi: float): self._set_fuel_ratio(phi=phi)
  def set_r_eq(self, r_eq: float): self._set_fuel_ratio(r_eq=r_eq)
  
  def set_pressure(self, pressure: float): self.pressure = pressure
  def set_density(self, density: float): self.density = density # for UV

  def set_materials(self, materials: List[Material]): self.materials = materials
  def set_massf(self, massf: bool): self.massf = massf
  def set_inserts(self, inserts: Union[str, List[Union[str, Material]]]): self.inserts = self._format_input_list(inserts)
  def set_omits(self, omits: Union[str, List[Union[str, Material]]]): self.omits = self._format_input_list(omits)

  def set_filename(self, filename: str):
    if ".inp" in filename or "/" in filename: # Must be a string of alphanumeric characters
      raise ValueError("Cannot save to filename with .inp or / in it")
    self.filename = filename
  
  def set_pressure_units(self, pressure_units: str):
    pressure_units = pressure_units.lower() # We only have a few options for pressure units
    if pressure_units not in OPTIONS_PRES_UNITS:
      raise ValueError("pressure unit must be in " + ", ".join(OPTIONS_PRES_UNITS))
    self.pressure_units = pressure_units
  
  # Density for UV problems
  def set_density_units(self, density_units: str):
    density_units = density_units.lower() # We only have a few options for density units
    if density_units not in OPTIONS_DENS_UNITS:
      raise ValueError("density unit must be in " + ", ".join(OPTIONS_DENS_UNITS))
    self.density_units = density_units
  
  def set_absolute_o_f(self) -> float:
    # Set an o_f ratio assuming that the wt_percent of each of our materials is an absolute percentage
    sum_ox   = sum([item.wt_percent for item in filter(lambda x: isinstance(x, Oxidizer), self.materials)])
    sum_fuel = sum([item.wt_percent for item in filter(lambda x: isinstance(x, Fuel), self.materials)])
    o_f = sum_ox/sum_fuel
    self.set_o_f(o_f)
    return o_f
  
  def run(self, *materials) -> Output:
    if self.ratio_name == None:
      raise TypeError("No reactant ratio specified, must set phi, or o/f, or %f, etc.")
    if len(materials) > 0: # If they specify materials, update our list
      self.materials = materials
    try:
      self.make_input_file(self.materials)
    except OSError: # Sometimes things like dropbox lock the file so we can't access it
      raise RuntimeError("unable to open input file for writing...")
    run_cea_backend(self.filename)
    return self.process_output()
  run_cea = run # alias for backward-compatibility
  
  def make_input_file(self, material_list: List[Material]): # chamber conditions and materials list
    # Make sure we have some materials
    if material_list is None or len(material_list) == 0:
      raise ValueError("must specify at least one Material in Problem")

    # Makes sure that all our materials are Material objects
    if any([not isinstance(x, Material) for x in material_list]):
      raise ValueError("all items in Problem material list must be Material objects")
    # First we need to check that all of our materials have either mole/wt fractions specified. No mixing.
    # If the first element has a wt_percent, then we check that all elements have a wt_percent
    if not all([(mat.is_mols() if material_list[0].is_mols() else mat.is_wt_percent()) for mat in material_list]):
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
  
    with open(self.filename+".inp", "w") as file:
      file.write("problem ")
      file.write(self.get_prefix_string())
      file.write("react  \n")
      # First we write all the fuels (I don't think there's a reason to write them in this order, but we'll do it anyway
      for fuel in fuels:
        file.write(fuel.get_CEA_str())
      # Then we write all the oxidizers
      for ox in oxidizers:
        file.write(ox.get_CEA_str())
      if self.inserts:
        file.write("insert "+" ".join(self.inserts)+"\n")
      if self.omits:
        file.write("omit "+" ".join(self.omits)+"\n")
        
      file.write("output {}trans\n".format("massf " if self.massf else "")) # Output is mass fraction is "massf"
      file.write(self.get_plt_string()) # output plotting string, if any
      file.write("end\n")
  
  def get_plt_string(self) -> str:
    return f"   plot {self.plt_keys:s}\n" # specify string formatting so None errors
  
  def _error_check_thermo_out_line(self, line):
    # Either raises an error or does nothing
    filePath = repr(self.filename+".out")
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
      if key in line:
        raise RuntimeError(prepend + value + filePath)
    
  
  def process_output(self) -> str:
    # Process Plot file
    out = Output()
    try:
      with open(self.filename+".plt", errors='ignore') as file:
        for line in file:
          new_line = line.split()
          if new_line[0] == '#':
            continue
          else:
            for key, value in zip(self.plt_keys.split(" "), new_line):
              out[key] = float(value)
    except FileNotFoundError:
      raise RuntimeError("CEA Failed to Run. Plot file wasn't generated for " + self.filename)
    return out

  ## THESE MUST BE IMPLEMENTED BY SUBCLASSES ##
  def get_prefix_string(self) -> str:
    raise NotImplementedError()


class DetonationProblem(Problem):
  problem_type = "det"
  plt_keys = "p t h mw cp gammas phi vel mach rho son cond pran"
  
  def get_prefix_string(self):
    toRet = []
    toRet.append("{}".format(self.problem_type))
    toRet.append("   p({}) = {:0.5f}".format(self.pressure_units, self.pressure))
    toRet.append("   {} = {:0.5f}".format(self.ratio_name, self.ratio_value))
    return "\n".join(toRet) + "\n"

  def process_output(self):
    out = Output()
    
    # We'll also open this file to get mass/mole fractions of all constituents and other properties
    with open(self.filename+".out") as file:
      # Ordered (in the file) list of terms that we're searching for
      search_terms = ["BURNED GAS", "DETONATION PARAMETERS", ("MASS FRACTIONS" if self.massf else "MOLE FRACTIONS")]
      search_i = 0
      in_search_area = False
      was_in_search_area = False # State variable so we increment search_i when it changes high to low
      passed_first_line = False # State variable so we ignore the empty line at the start of each section
      
      # The format of each section is "HEADING\n [empty line]\n [lines of data...]\n [empty line]
      #   So when we get a heading, we ignore the first line after, and keep going until the code for that section says to stop
      
      out.prod_c = Output()
      
      for line in file:
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
            split = re.findall("\S+", line) # match by sections which are not whitespace
            key = split[0].lstrip("*") # remove any asterisks
            out.prod_c[key] = float(split[1])
          ### Detonation Products ###
          elif search_i==1:
            if not line.strip():
              in_search_area = False
              continue
            split = re.findall("\S+", line) # match by sections which are not whitespace
            key = split[0].replace("/","_").lower() # convert to usable format
            if key == "det": # This means we've gone past the p, t, m, rho / initial lines
              in_search_area = False
              continue # break out
            out[key] = float(split[1])
          ### Output gas products ###
          elif search_i == 0:
            if not line.strip(): # If this line is empty, skip it
              continue
            split = re.findall("\S+", line) # match by sections which are not whitespace
            key = split[0]
            if key == "(dLV/dLP)t":
              out.dLV_dLP_t = float(split[1])
            elif key == "(dLV/dLT)p":
              out.dLV_dLT_p = float(split[1])
            elif key == "VISC,MILLIPOISE": # Keep going until we find the viscosity line
              out.visc = float(split[1]) * 0.0001 # convert millipoise to Pa-s
              in_search_area = False
              continue
    
    outPlt = super().process_output() # Process the plt file second so we can read errors listed in output file
    out.update(outPlt)
    
    # Also include the actual gamma and not just the isentropic gamma
    out.gamma = out.gammas*-out.dLV_dLP_t
    
    return out
    
class HPProblem(Problem):
  problem_type = "hp"
  plt_keys = "p t h mw cp gammas phi rho son cond pran"
  def get_prefix_string(self):
    toRet = []
    toRet.append("{}".format(self.problem_type))
    toRet.append("   p({}) = {:0.5f}".format(self.pressure_units, self.pressure))
    toRet.append("   {} = {:0.5f}".format(self.ratio_name, self.ratio_value))
    return "\n".join(toRet) + "\n"
  
  def process_output(self):
    out = Output()
    
    # We'll also open this file to get mass/mole fractions of all constituents and other properties
    with open(self.filename+".out") as file:
      # Ordered (in the file) list of terms that we're searching for
      search_terms = ["THERMODYNAMIC PROPERTIES", ("MASS FRACTIONS" if self.massf else "MOLE FRACTIONS")]
      search_i = 0
      in_search_area = False
      was_in_search_area = False # State variable so we increment search_i when it changes high to low
      passed_first_line = False # State variable so we ignore the empty line at the start of each section
      
      # The format of each section is "HEADING\n [empty line]\n [lines of data...]\n [empty line]
      #   So when we get a heading, we ignore the first line after, and keep going until the code for that section says to stop
      
      out.prod_c = Output()
      
      for line in file:
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
            split = re.findall("\S+", line) # match by sections which are not whitespace
            key = split[0].lstrip("*") # remove any asterisks
            out.prod_c[key] = float(split[1])
          ### Output gas products ###
          elif search_i == 0:
            if not line.strip(): # If this line is empty, skip it
              continue
            split = re.findall("\S+", line) # match by sections which are not whitespace
            key = split[0]
            if key == "(dLV/dLP)t":
              out.dLV_dLP_t = float(split[1])
            elif key == "(dLV/dLT)p":
              out.dLV_dLT_p = float(split[1])
            elif key == "VISC,MILLIPOISE": # Keep going until we find the viscosity line
              out.visc = float(split[1]) * 0.0001 # convert millipoise to Pa-s
              in_search_area = False
              continue
    
    outPlt = super().process_output() # Process the plt file second so we can read errors listed in output file
    out.update(outPlt)
    
    # Also include the actual gamma and not just the isentropic gamma
    out.gamma = out.gammas*-out.dLV_dLP_t
    
    return out

class TPMaterial(Material):
  def __init__(self, parent):
    try:
      p = parent
      self.parent = p
      self.output_type = self.parent.output_type
      self.is_fuel = self.parent.is_fuel
      super().__init__(p.name, p.temp, p.wt_percent, p.mols, p.chemical_composition, p.hf)
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
    self.temperature_units = None
    self.set_temperature_units(temperature_units)
  
  def set_temperature(self, temperature: float): self.temperature = temperature
  
  def set_temperature_units(self, temperature_units: str):
    temperature_units = temperature_units.lower()  # We only have a few options for pressure units
    if temperature_units not in OPTIONS_TEMP_UNITS:
      raise ValueError("pressure unit must be in " + ", ".join(OPTIONS_TEMP_UNITS))
    self.temperature_units = temperature_units

  def get_prefix_string(self):
    toRet = []
    toRet.append("{}".format(self.problem_type))
    toRet.append("   p({}) = {:0.5f}".format(self.pressure_units, self.pressure))
    toRet.append("   t({}) = {:0.5f}".format(self.temperature_units, self.temperature))
    toRet.append("   {} = {:0.5f}".format(self.ratio_name, self.ratio_value))
    return "\n".join(toRet) + "\n"
  
  def make_input_file(self, material_list: List[Material]): # chamber conditions and materials list
    # Here we convert all the materials to TPMaterials before calling the main function, so that each one doesn't report temperature
    new_material_list = [TPMaterial(material) for material in material_list]
    return super().make_input_file(new_material_list)

class RocketProblem(Problem):
  problem_type = "rocket"
  plt_keys = "p t isp ivac m mw cp gam o/f cf rho son mach phi h cond pran ae pip"
    
  def __init__(self, *args, sup: float=None, sub: float=None, ae_at: float=None, pip: float=None, analysis_type:str="equilibrium", fac_ac:float=None, fac_ma:float=None, nfz: int=None, custom_nfz: float=None,
               **kwargs):
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
    self.nozzle_ratio_value = [pip if pip else (sup if sup else sub), ] # 1-element array
    self.fac_type = None
    self.fac_value = None
    self.analysis_type = None  # Specify before calling
    if fac_ac: self.set_fac_ac(fac_ac)
    if fac_ma: self.set_fac_ma(fac_ma)
    self.set_analysis_type(analysis_type, nfz, custom_nfz=custom_nfz)
    
  
  
  def set_sup(self, sup): self.nozzle_ratio_name = "sup"; self.nozzle_ratio_value[-1] = sup
  def set_sub(self, sub): self.nozzle_ratio_name = "sub"; self.nozzle_ratio_value[-1] = sub
  def set_ae_at(self, ae_at): self.set_sup(ae_at)
  
  def set_pip(self, pip): self.nozzle_ratio_name = "pip"; self.nozzle_ratio_value[-1] = pip
  
  def set_analysis_type(self, analysis_type, nfz=None, custom_nfz=None):
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
      toPut = type(self).plt_keys # Get copy from class instance
      for orig in ("isp", "ivac", "cf", "mach"):
        toPut = toPut.replace(" "+orig, " "+orig+"fz", 1) # Add 'fz' to each of these properties
      self.plt_keys = toPut
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
    toRet = list()
    toRet.append("{} {}".format(self.problem_type, self.analysis_type))
    toRet.append("   p({}) = {:0.5f}".format(self.pressure_units, self.pressure))
    toRet.append("   {} = {:0.5f}".format(self.ratio_name, self.ratio_value))
    # enable using multiple area ratio for custom frozen point
    toRet.append("   {} = ".format(self.nozzle_ratio_name) + ",".join("{:0.5f}".format(rat) for rat in self.nozzle_ratio_value))
    if self.fac_type:
      toRet.append("   fac {} = {:0.5f}".format(self.fac_type, self.fac_value))
    return "\n".join(toRet) + "\n"
  
  def process_output(self):
    out = Output()
    
    # We'll also open this file to get mass/mole fractions of all constituents and other properties
    with open(self.filename+".out") as file:
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
      def flMap(splitline, key):
        if has_fac: # If finite area combustor, we have more columns
          mapping = {"c": 1, "f": 2, "t": 3, "e": end_col}
        else: # Otherwise, normal
          mapping = {"c": 1, "t": 2, "e": end_col}
        return float(splitline[mapping[key]])
      
      for line in file:
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
            split = re.findall("\S+", line) # match by sections which are not whitespace
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
            split = re.findall("\S+", line) # match by sections which are not whitespace
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
            split = re.findall("\S+", line) # match by sections which are not whitespace
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
    
    try:
      chamber_count = 1
      if has_fac:
        fac_count = 2
        throat_count = 3
        exit_count = 4
      else:
        fac_count = -1
        throat_count = 2
        exit_count = end_col
      
      with open(self.filename+".plt", errors='ignore') as file:
        for counter, line in enumerate(file):
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
            
    except FileNotFoundError:
      raise RuntimeError("CEA Failed to Run. Plot file wasn't generated for " + self.filename)
            

    if "frozen" in self.analysis_type:
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

    return out

#// Added for UV Implementation
#// Mostly just a copy-paste with some changes
#//  moving it from pressure-defined to density-defined
class UVProblem(Problem):
  problem_type = "uv"
  plt_keys = "p t rho h u g s m mw cp gammas phi rho son cond pran"
  
  def __init__(self, *args, density: float=1, # Chamber reactants relative volume
               density_units: str="kg", **kwargs):
    self.density = density
    self.density_units = density_units
    self.set_density_units(density_units)
    super().__init__(*args, **kwargs)

  def get_prefix_string(self):
    toRet = []
    toRet.append("{}".format(self.problem_type))
    toRet.append("   rho({}) = {:0.5f}".format(self.density_units, self.density))
    toRet.append("   {} = {:0.5f}".format(self.ratio_name, self.ratio_value))
    return "\n".join(toRet) + "\n"
  
  def process_output(self):
    out = Output()
    
    # We'll also open this file to get mass/mole fractions of all constituents and other properties
    with open(self.filename+".out") as file:
      # Ordered (in the file) list of terms that we're searching for
      search_terms = ["THERMODYNAMIC PROPERTIES", ("MASS FRACTIONS" if self.massf else "MOLE FRACTIONS")]
      search_i = 0
      in_search_area = False
      was_in_search_area = False # State variable so we increment search_i when it changes high to low
      passed_first_line = False # State variable so we ignore the empty line at the start of each section
      
      # The format of each section is "HEADING\n [empty line]\n [lines of data...]\n [empty line]
      #   So when we get a heading, we ignore the first line after, and keep going until the code for that section says to stop
      
      out.prod_c = Output()
      
      for line in file:
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
            split = re.findall("\S+", line) # match by sections which are not whitespace
            key = split[0].lstrip("*") # remove any asterisks
            out.prod_c[key] = float(split[1])
          ### Output gas products ###
          elif search_i == 0:
            if not line.strip(): # If this line is empty, skip it
              continue
            split = re.findall("\S+", line) # match by sections which are not whitespace
            key = split[0]
            if key == "(dLV/dLP)t":
              out.dLV_dLP_t = float(split[1])
            elif key == "(dLV/dLT)p":
              out.dLV_dLT_p = float(split[1])
            elif key == "VISC,MILLIPOISE": # Keep going until we find the viscosity line
              out.visc = float(split[1]) * 0.0001 # convert millipoise to Pa-s
              in_search_area = False
              continue
    
    outPlt = super().process_output() # Process the plt file second so we can read errors listed in output file
    out.update(outPlt)
    
    # Also include the actual gamma and not just the isentropic gamma
    out.gamma = out.gammas*-out.dLV_dLP_t
    
    return out


if __name__ == "__main__":
  # Example: Set up a simple Water + Peroxide reaction
  water = Fuel("H2O(L)")
  hPeroxide = Oxidizer("H2O2(L)")
  
  chamber_pressure = 1000
  exit_pressure = 14.7
  
  problem = RocketProblem(pressure=1000, o_f=5, omits="Al(cr)")
  
  for ae_at in [chamber_pressure/exit_pressure, 80]:
    problem.set_ae_at(ae_at)
    data = problem.run_cea(water, hPeroxide)
    print(data)
