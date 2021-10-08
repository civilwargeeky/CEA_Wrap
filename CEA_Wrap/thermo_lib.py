import re, os.path
from time import time
from .utils import _get_asset, Output

class ThermoMaterial(Output):
  def __init__(self, *args, **kwargs):
    super().__init__()
    
    required = (
                "name", # Name of the material
                "reference", # Reference, if given. Otherwise ""
                "elements", # Dictionary of element: numerical value specified
                "condensed", # True if condensed phase, False otherwise
                "mol_wt", # Molecular weight in g/mol
                "hf", # float heat of formation at 298.15 in kJ/mol (or assigned enthalpy if 0 temp range)
                "temp_ranges", # List of two-tuples of [range start, range end] (K)
                "reactant_only", # True if material only shows up in reactants, False otherwise
               )
    
    for elem in required:
      try:
        self[elem] = kwargs[elem]
      except KeyError:
        raise KeyError("ThermoMaterial did not receive keyword arg '{}'".format(elem)) from None
        
  def defined_at(temp):
    pass

def load_thermo_file(filename = _get_asset("thermo_spg.inp")):
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
      mol_wt = mw,
      hf = hf,
      temp_ranges = temps,
      reactant_only = after_air,
    )
    
  materials = {}
  after_air = False # All materials after and including Air may only be defined as reactants (they don't show up in products)
  # Source: The CEA specification
  
  print("Processing {} thermo input file".format(os.path.basename(filename)))
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
            materials[mat.name] = mat
          except AssertionError: # errors on "thermo" line because it isn't a proper line
            pass
          cur_mat_lines.clear() # clear this array for further use
        if line == "END REACTANTS": # Specified at the end of file
          break
      # Regardless of whether we are at definition line or not
      cur_mat_lines.append(line)
  elapsed = time() - start_time
  print("Processed {} materials in {:0.4f} s ==> {:0.2f} materials/s".format(len(materials), elapsed, len(materials)/elapsed))
  return materials
  
if __name__ == "__main__":
  load_thermo_file()