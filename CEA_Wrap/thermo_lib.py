from dataclasses import dataclass
import difflib  # for nearby element matches
import logging
import os.path
from time import time

from .mydifflib import get_close_matches_indexes as close_match_index  # for nearby element matches
from .utils import Output, get_asset, get_thermo_inp_location

log = logging.getLogger(__name__)

# These are molecules in the default thermo.inp that have multiple entries specified
# If entries besides these are defined multiple times, we should give user a warning
DEFAULT_DOUBLED_NAMES = ["Co(b)", "Cr(cr)", "Cr2O3(I)", "Cr2O3(I)", "Fe(a)", "Fe2O3(cr)", "Fe3O4(cr)", "K2S(cr)", "Na2S(cr)", "Ni(cr)", "SnS(cr)"]

@dataclass(frozen=True)
class ThermoMaterial:
  """
  A ThermoMaterial represents a single entry in a thermo.inp file

  :param name: Name of the material
  :param reference: Reference, if given. Otherwise ""
  :param elements: Dictionary of element: numerical value specified
  :param condensed: True if condensed phase, False otherwise
  :param mw: Molecular weight in g/mol
  :param hf: float heat of formation at 298.15 in kJ/mol (or assigned enthalpy if 0 temp range)
  :param temp_ranges: List of two-tuples of [range start, range end] (K)
  :param reactant_only: True if material only shows up in reactants, False otherwise
  """
  name: str
  reference: str
  elements: dict[str, int]
  condensed: bool
  mw: float
  hf: float
  temp_ranges: list[tuple[float, float]]
  reactant_only: bool

  def __getitem__(self, name:str): # Allow this to act as a dictionary also
    return getattr(self, name)

  def defined_at(self, temp:float) -> bool:
    # Returns True if the element is defined in any of the reactant's temperature ranges
    for low, high in self.temp_ranges:
      if low <= temp <= high:
        return True
    return False

class ThermoInterface:
  def __init__(self, filename:None|str=None):
    """
    Instantiates a ThermoInterface that parses a given thermo.inp file (note: does not need to be called "thermo.inp")

    :param filename: The path to the thermo.inp file. If None, uses the default path according to get_asset, defaults to None
    """
    self.filename = None
    self.thermo_materials = None
    self.load(get_thermo_inp_location() if filename is None else filename)
  
  def load(self, filename:str):
    self.filename = filename
    self.thermo_materials = self.read_thermo_file(filename)
  
  def read_thermo_file(self, filename:str) -> dict[str, ThermoMaterial]:
    """
      Loads the thermo input file at filename
      Returns a dictionary of name: ThermoMaterial representing the material
    """
    split_by_num = lambda string, num: [string[i:i+num] for i in range(0, len(string), num)]
    
    def process_lines(lines, after_air):
      # Function to process a list of lines and return a ThermoMaterial
      # Process First Line: Name and reference
      name = lines[0].split(" ")[0] # first word of first line
      if name == "thermo":
        raise AssertionError("can't process 'thermo' line")
      ref = " ".join(lines[0].split()[1:])
      
      # Process second line: number of temperature ranges, elements, mw, heat of formation, condensed or not
      num_ranges = int(lines[1].split()[0]) # number of temperature ranges
      elems_string = lines[1][10:50] # these columns for elements
      # each element gets 8 characters, 2 letters and 6.2F
      # split by 8 characters
      elems_strs = split_by_num(elems_string, 8)
      # ... what happens if something has more than 5 elements?
      elems = {}
      for elem in elems_strs:
        if elem[0].isspace(): # if first is whitespace, no more elements
          break
        elems[elem[0:2].rstrip()] = float(elem[2:]) # dict of element and float amount
      condensed = lines[1][51] != "0" # "Zero for gas and non-zero for condensed phase"
      mw = float(lines[1][52:65])
      hf = float(lines[1][65:80])
      
      # Process Temperature Ranges
      if num_ranges == 0: # Only specified at one temperature +- 10K
        base_temp = float(lines[2][1:11])
        temps = [(base_temp-10, base_temp+10)]
      else:
        temps = []
        for line in lines[2::3]: # get every third line from 2 to end. Contains temp range
          low_temp  = float(line[1:11])
          high_temp = float(line[11:21])
          temps.append((low_temp, high_temp))
          
      return ThermoMaterial(
        name = name, 
        reference = ref,
        elements = elems,
        condensed = condensed,
        mw = mw,
        hf = hf,
        temp_ranges = temps,
        reactant_only = after_air,
      )
      
    materials = {}
    after_air = False # All materials after and including Air may only be defined as reactants (they don't show up in products)
    # Source: The CEA specification
    
    log.debug("Processing {} thermo input file".format(os.path.basename(filename)))
    start_time = time()
    
    with open(filename) as file:
      cur_mat_lines = []
      for line in file:
        line = line.rstrip() # only rstrip so I can index lines
        if line[0] == "!":
          continue
        if line[0].isalpha() or line[0] in ["(",]: # If starts with a letter -> a material == a line for specifying a new material
          if line == "END PRODUCTS": # all after this will be reactant only
            after_air = True
            continue
          if len(cur_mat_lines) > 0: # We will now process the previous material that we are at a new one
            try:
              mat = process_lines(cur_mat_lines, after_air)
              if mat.name in materials: # For some God-forsaken reason, some materials have multiple entries with identical properties for multiple temp ranges
                log.debug("Name overlap for Material '{}'".format(mat.name))
                if mat.name not in DEFAULT_DOUBLED_NAMES:
                  log.warning(f"User material '{mat.name}' shares name with default material. Old material not overwritten. Material temperature ranges appended") 
                materials[mat.name].temp_ranges.extend(mat.temp_ranges)
              else: # First time we're seeing it
                materials[mat.name] = mat
            except AssertionError: # errors on "thermo" line because it isn't a proper line
              pass
            cur_mat_lines.clear() # clear this array for further use
          if line == "END REACTANTS": # Specified at the end of file
            break
        # Regardless of whether we are at definition line or not
        cur_mat_lines.append(line)
    elapsed = time() - start_time
    log.debug("Processed {} materials in {:0.4f} s ==> {:0.2f} materials/s".format(len(materials), elapsed, len(materials)/elapsed))
    return materials
  
  def get_close_matches(self, name:str, n:int=6) -> list[str]:
    """
    Returns names in the thermo library that are close matches for the given name (using difflib)

    :param name: The name to search for close matches to
    :param n: Maximum number of close matches, defaults to 6
    :return: Between 0 and `n` close matches for the given name within the names in the thermo library
    """
    n1 = int(n/2) # round down
    n2 = int(round(n/2,0)) # round up
    name_list = list(self.thermo_materials) # maintain ordering for indexing
    # Get matches for the material itself. Matches things like "CH5" --> "CH4" or "Al(cr)" --> "AL(cr)"
    close_matches1 = difflib.get_close_matches(name, name_list, n=n1)
    # Then get matches for partially filled names. Matches things like "C8H18" --> "C8H18,isohexalate" or whatever
    close_matches2_ind = close_match_index(name, list(map(lambda x: x[:len(name)], name_list)), n=n2)
    close_matches2 = [name_list[i] for i in close_matches2_ind]
    return list(sorted(set(close_matches1 + close_matches2))) # set so we have unique entries, then sort and convert to list
  
  def __getattr__(self, name):
    """ Allows for item access and all dict methods of thermo materials dict """
    return getattr(self.thermo_materials, name)


def get_simple_thermo_lib_line(name: str, comment: str, atoms: dict[str: int], is_condensed: bool, 
                                 mol_wt: float, hf: float, temperature: float=298.150) -> str:
  """
  Returns and prints a string which can be inserted to represent a molecule in the ThermoLib
  Any entries which are longer than the allotted space results in an error

  :param name: 18 chars max, Species name, such as CO or Fe2O3. These are assumed to be in a gas phase
    unless appended with (a) for aqueous solution, (cr) for crystalline (solid), or (L) for liquid phase.
  :param comment: 62 char max, Human-readable name and additional information
  :param atoms: 5 atom max, A dictionary of "atom symbol": number of atoms in molecule.
    Format: {"atomic symbol":number of atoms, "atomic symbol" : number of atoms, etc.}
    Example: H2O would be {"H": 2, "O": 1}
    Special value is "E" which represents an electron for ionic compounds
  :param is_condensed: True if condensed phase (non-gas), False otherwise
  :param mol_wt: Molecular weight of molecule in g/mol or kg/kmol
  :param hf: Heat of formation at 298.15 K, in J/mol
  :param temperature: Default 298.150, the temperature that these heats of formation apply to.
  """
  if len(name) > 18:
    raise ValueError("Name field is {} characters long, which is more than 18".format(len(name)))
  if " " in name:
    raise ValueError("Name field cannot contain any spaces")
  if len(comment) > 62:
    raise ValueError("Comment field is {} characters long, which is more than 62".format(len(name)))
  if len(atoms) > 5:
    raise ValueError("Got {} atoms in molecule. Can only support 5".format(len(atoms)))
  atomString = "    0.00"*5 # Initialize with all the empty atom portions
  for atom_name, atom_num in atoms.items():
    if len(atom_name) > 2:
      raise ValueError("Atom names must be 2 chars or less")
    if atom_num > 999:
      raise ValueError("Number of {} atoms must be less than 999".format(atom_name))
    if atom_num <= 0 and atom_name != "E": # Special value: Positively charged ions have "E: -1"
      raise ValueError("Number of {} atoms must be greater than 0".format(atom_name))
    atomString = "{:<2}{:6.2f}".format(atom_name.upper(), atom_num) + atomString
  atomString = atomString[:40] # Then truncate to 40 chars to remove extra empty atom portions
  if len("{:13.5f}".format(mol_wt)) > 13:
    raise ValueError("molecular weight has too many digits")
  if len("{:15.5f}".format(hf)) > 15:
    raise ValueError("heat of formation has too many digits")
  
  string = "{:<18}{:<62}\r\n 0        {:40} {:1}{: 13.5f}{: 15.5f}\r\n".format(
    name,
    comment,
    atomString,
    "1" if is_condensed else "0",
    mol_wt,
    hf,
  )
  # Add in the temperature line that everyone forgets
  string+="\r\n {:10.3f}      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000".format(temperature)
  return string

def print_simple_thermo_lib_line(*args, **kwargs):
  print(get_simple_thermo_lib_line(*args, **kwargs))