import subprocess, re, os.path, shutil
import importlib.resources

with importlib.resources.path(__package__, "FCEA2.exe") as manager:
  CEA_LOCATION = str(manager)

for file in ["thermo.lib", "trans.lib"]:
  if not os.path.isfile(file):
    print(file+" not found in current directory. Copying from package...")
    with importlib.resources.path(__package__, file) as manager:
      shutil.copyfile(manager, file)

class Material:
  output_type = None # MUST BE SPECIFIED IN SUBCLASSES
  
  def __init__(self, name, temp=298, wt_percent=None, mols=None, chemical_composition = None, hf = None):
    if wt_percent is None and mols is None: # If neither is specified, user probably doesn't care
      wt_percent = 100
    
    self.name = name # Name to input into CEA
    self.temp = temp # temperature in K, defaults to 298K
    self.wt_percent = wt_percent # Should be set when actually running through
    self.mols = mols # Mole ratio for this material
    self.chemical_composition = chemical_composition # chemical composition such as "LI 1 B 1 H 4" for LiBH4. If defined, will not use CEA default values
    self.hf = hf # Enthalpy of formation, if needs to be defined.
    if chemical_composition != None and hf == None:
      raise ValueError("Elements entered molecule by molecule must have defined hf")
    
    if wt_percent and mols:
      raise TypeError("Material cannot have both wt_percent and mols specified")
  
  def set_wt_percent(self, wt_percent):
    self.mols = None # Can only have one
    self.wt_percent = wt_percent
  
  def set_mols(self, mols):
    self.wt_percent = None # Can only have one
    self.mols = mols
    
  def is_mols(self): # Helper function for a material being in wt_percent or mols
    return self.mols is not None
    
  def is_wt_percent(self): # Helper function for a material being in wt_percent or mols
    return self.wt_percent is not None
    
  def get_CEA_str(self):
    # Specify whether to use the str/val for weight or mols.
    name, ratio = "wt" if self.wt_percent is not None else "mol", self.wt_percent if self.wt_percent is not None else self.mols
    if ratio < 0:
      raise ValueError("Cannot have {} with <0 {} percent!".format(self.name, name))
    elif ratio == 0: # This allows us to not include materials that are set to 0 percent
      return ""
    else:
      string = "   {}={} {}={:0.5f}  t,k={:0.5f}".format(self.output_type, self.name, name, ratio, self.temp)
      if self.hf:
        string += " h,kj/mol={:0.5}".format(self.hf)
      if self.chemical_composition:
        string += " " + self.chemical_composition
      return string + " \n"
      
class Fuel(Material):
  output_type = "fuel"
F = Fuel # Alias

class Oxidizer(Material):
  output_type = "ox"
O = Oxidizer # Alias

## Plan: Can also have composite Fuel/Oxidizer made up of percentages of other components

def run_cea_backend(filename):
  ret = subprocess.run(CEA_LOCATION, input=filename+"\n", text=True, stdout=subprocess.DEVNULL)
  if ret.returncode != 0:
    print(ret)
    raise RuntimeError("Running CEA failed with errors")
  return ret
  
class Output(dict): # This is just a dictionary that you can also use dot notation to access
  def __init__(self): # Explicitly must receive no arguments, because I don't want to deal with constructor properties
    super()
  def __getattr__(self, name):
    return self[name]
  def __setattr__(self, name, value):
    if name.startswith("_"):
      super().__setattr__(name, value)
    else:
      self[name] = value

# My idea is that we could have different types of problems with similar methods for like get_prefix_string and things
class Problem:
  ### NOTE : ALL PROBLEM SUBCLASSES MUST SPECIFY PROBLEM TYPE AND PLOT KEYS ###
  problem_type = None
  plt_keys = None # plt_keys should be a space-separated string of items to go into the "plt" command of FCEA2
  
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
        
  
  # All arguments must be specified by keyword
  def __init__(self, *, pressure=1000, materials=None, massf=False, filename="my_output", pressure_units="psi", inserts=None, omits=None, **kwargs):
    self.massf = massf
    self.materials = materials
    self.inserts = inserts
    self.omits = omits
    self.pressure = pressure
    self.set_pressure_units(pressure_units)
    self.set_filename(filename)
    
    self.ratio_name = None
    self.ratio_value = None
    self._set_fuel_ratio(**kwargs)
    
    diff = set(kwargs).difference(set(self._ratio_options))
    if diff: # If there are any kwargs keys that aren't in _ratio_options keys, we should error
      raise TypeError(self.__class__.__name__+"() got an unexpected keyword argument(s): " + ",".join(diff))
    
  def set_p_f(self, p_f): self._set_fuel_ratio(p_f=p_f)
  def set_f_o(self, f_o): self._set_fuel_ratio(f_o=f_o)
  def set_o_f(self, o_f): self._set_fuel_ratio(o_f=o_f)
  def set_phi(self, phi): self._set_fuel_ratio(phi=phi)
  def set_r_eq(self, r_eq): self._set_fuel_ratio(r_eq=r_eq)
  
  def set_pressure(self, pressure): self.pressure = pressure
  def set_materials(self, materials): self.materials = materials
  def set_massf(self, massf): self.massf = massf
  def set_filename(self, filename): 
    if ".inp" in filename or "/" in filename: # Must be a string of alphanumeric characters
      raise ValueError("Cannot save to filename with .inp or / in it")
    self.filename = filename
  def set_pressure_units(self, pressure_units):
    pressure_units = pressure_units.lower() # We only have a few options for pressure units
    choices = ["bar", "atm", "psi", "mmh"]
    if pressure_units not in choices:
      raise ValueError("pressure unit must be in " + ", ".join(choices))
    self.pressure_units = pressure_units
  
  def set_absolute_o_f(self):
    # Set an o_f ratio assuming that the wt_percent of each of our materials is an absolute percentage
    sum_ox   = sum([item.wt_percent for item in filter(lambda x: isinstance(x, Oxidizer), self.materials)])
    sum_fuel = sum([item.wt_percent for item in filter(lambda x: isinstance(x, Fuel), self.materials)])
    o_f = sum_ox/sum_fuel
    self.set_o_f(o_f)
    return o_f
  
  def run_cea(self, *materials):
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
  
  def make_input_file(self, material_list): # chamber conditions and materials list
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
      
    fuels = list(filter(lambda x: isinstance(x, Fuel), material_list))
    oxidizers = list(filter(lambda x: isinstance(x, Oxidizer), material_list))
    
    # Here we check for monopropellant cases. In a monopropellant case, you must specify only fuels and an o_f of 0 or you get weird results
    if len(fuels) == 0 and len(oxidizers) > 0:
      raise ValueError("Monopropellant problems must only specify fuels, not oxidizers")
    if len(oxidizers) == 0 and (self.ratio_name != "o/f" or self.ratio_value != 0):
      raise ValueError("Monopropellant problems must be run with o/f of 0")
  
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
  
  def get_plt_string(self):
    return f"   plot {self.plt_keys:s}\n" # specify string formatting so None errors
  
  def process_output(self):
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
  def get_prefix_string(self):
    raise NotImplementedError()

class Detonation_Problem(Problem):
  problem_type = "det"
  plt_keys = "p t h mw cp gammas phi vel mach rho son"
  
  def get_prefix_string(self):
    toRet = []
    toRet.append("{}".format(self.problem_type))
    toRet.append("   p({}) = {:0.5f}".format(self.pressure_units, self.pressure))
    toRet.append("   {} = {:0.5f}".format(self.ratio_name, self.ratio_value))
    return "\n".join(toRet) + "\n"

  def process_output(self):
    out = super().process_output() # Process the plt file first

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
        if "FATAL" in line:
          raise RuntimeError("CEA Failed to Run. FATAL error in input/output file: " + self.filename)
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
    
    # Also include the actual gamma and not just the isentropic gamma
    out.gamma = out.gammas*-out.dLV_dLP_t
    
    return out
    
class HP_Problem(Problem):
  problem_type = "hp"
  plt_keys = "p t h mw cp gammas phi rho son"
  def get_prefix_string(self):
    toRet = []
    toRet.append("{}".format(self.problem_type))
    toRet.append("   p({}) = {:0.5f}".format(self.pressure_units, self.pressure))
    toRet.append("   {} = {:0.5f}".format(self.ratio_name, self.ratio_value))
    return "\n".join(toRet) + "\n"
  
  def process_output(self):
    out = super().process_output() # Process the plt file first
  
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
        if "FATAL" in line:
          raise RuntimeError("CEA Failed to Run. FATAL error in input/output file: " + self.filename)
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
    
    # Also include the actual gamma and not just the isentropic gamma
    out.gamma = out.gammas*-out.dLV_dLP_t
    
    return out

class Rocket_Problem(Problem):
  problem_type = "rocket"
  plt_keys = "p t isp ivac m mw cp gam o/f cf rho son mach phi h"
    
  def __init__(self, *args, sup=None, sub=None, ae_at=None, analysis_type="equilibrium", **kwargs):
    super().__init__(*args, **kwargs)
    
    if ae_at:
      sup = ae_at
    if not sup and not sub: # Default if nothing is specified
      sup = 1
    if sup and sub:
      raise ValueError("Can only specify supersonic or subsonic area ratio, not both")

    analysis_type = analysis_type.lower() # ensure case because we check for frozen by literal
    self.area_ratio_name = "sup" if sup else "sub" # Can only specify one supersonic or subsonic area ratio
    self.area_ratio_value = sup if sup else sub
    self.ae_at = ae_at
    self.analysis_type = analysis_type # equilibrium or frozen
    
    if "equilibrium" in analysis_type and "frozen" in analysis_type:
      raise ValueError("Rocket_Problem does not support combined equilibrium-frozen calculations")
  
  
  def set_sup(self, sup): self.area_ratio_name = "sup"; self.area_ratio_value = sup
  def set_sub(self, sub): self.area_ratio_name = "sub"; self.area_ratio_value = sub
  def set_ae_at(self, ae_at): self.set_sup(ae_at)
  
  def get_prefix_string(self):
    toRet = []
    toRet.append("{} {}".format(self.problem_type, self.analysis_type))
    toRet.append("   p({}) = {:0.5f}".format(self.pressure_units, self.pressure))
    toRet.append("   {} = {:0.5f}".format(self.ratio_name, self.ratio_value))
    toRet.append("   {} = {:0.5f}".format(self.area_ratio_name, self.area_ratio_value)) # Actually wait this is dumb.... should be sup and sub separate
    return "\n".join(toRet) + "\n"
  
  def process_output(self):
    out = Output()
    try:
      with open(self.filename+".plt", errors='ignore') as file:
        for counter, line in enumerate(file):
          new_line = line.split()
          if counter == 0:
            continue
          elif  counter == 1: # chamber properties
            # Values common to all
            out.o_f = float(new_line[8])
            out.phi = float(new_line[13])
            
            # Section Properties 
            out.c_p = float(new_line[0])
            out.c_t = float(new_line[1])
            out.c_m = float(new_line[4])
            out.c_mw = float(new_line[5])
            out.c_cp = float(new_line[6])
            out.c_gammas = float(new_line[7])
            out.c_rho = float(new_line[10])
            out.c_son = float(new_line[11])
            out.c_h = float(new_line[14])
          elif counter == 2: # throat properties
            # Throat/Nozzle Values
            out.t_isp = float(new_line[2])/9.81
            out.t_ivac = float(new_line[3])/9.81
            out.t_cf = float(new_line[9])
            
            # Section Properties  
            out.t_p = float(new_line[0])
            out.t_t = float(new_line[1])
            out.t_m = float(new_line[4])
            out.t_mw = float(new_line[5])
            out.t_cp = float(new_line[6])
            out.t_gammas = float(new_line[7])
            out.t_rho = float(new_line[10])
            out.t_son = float(new_line[11])
            out.t_h = float(new_line[14])
          elif counter == 3: # noz properties
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
            out.cp = float(new_line[6])
            out.gammas = float(new_line[7])
            out.rho = float(new_line[10])
            out.son = float(new_line[11])
            out.h = float(new_line[14])
    except FileNotFoundError:
      raise RuntimeError("CEA Failed to Run. Plot file wasn't generated for " + self.filename)
            
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
      
      out.prod_c = Output()
      out.prod_t = Output()
      out.prod_e = Output()
      
      for line in file:
        if "FATAL" in line:
          raise RuntimeError("CEA Failed to Run. FATAL error in input/output file: " + self.filename)
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
              out.prod_c[key] = float(split[1])
              out.prod_t[key] = float(split[2])
              out.prod_e[key] = float(split[3])
          ### Output gas products ###
          elif search_i == 0:
            if not line.strip(): # If this line is empty, skip it
              continue
            split = re.findall("\S+", line) # match by sections which are not whitespace
            key = split[0]
            if key == "(dLV/dLP)t": # Will not exist in frozen
              out.c_dLV_dLP_t = float(split[1])
              out.t_dLV_dLP_t = float(split[2])
              out.dLV_dLP_t = float(split[3])
            elif key == "(dLV/dLT)p":# Will not exist in frozen
              out.c_dLV_dLT_p = float(split[1])
              out.t_dLV_dLT_p = float(split[2])
              out.dLV_dLT_p = float(split[3])
            elif key == "VISC,MILLIPOISE": # Keep going until we find the viscosity line
              out.c_visc = float(split[1]) * 0.0001 # convert millipoise to Pa-s
              out.t_visc = float(split[2]) * 0.0001 # convert millipoise to Pa-s
              out.visc = float(split[3]) * 0.0001 # convert millipoise to Pa-s
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

class Data_Collector(Output):
  def __init__(self, *args, keys=[], chamber_keys=[], exit_keys=[]):
    if len(args) > 0:
      if keys:
        raise TypeError("Data_Collector should not receive both a list of arguments and the 'keys' keyword")
      keys = list(args)
    self._prototype = list
    self._add_element = list.append
    self._keys = keys
    self._chamber_keys = chamber_keys
    self._exit_keys = exit_keys
    
    if set(chamber_keys).intersection(exit_keys):
      raise ValueError("Can't have product keys in both chamber and exit")
      
    if not all([isinstance(val, str) for val in self._keys + self._chamber_keys + self._exit_keys]):
      raise ValueError("All Data_Collector keys must be strings")
      
    for key in keys + chamber_keys + exit_keys:
      self[key] = self._prototype([])
    
  def add_data(self, data):
    for key in self._keys:
      self._add_element(self[key], data[key])
    for key in self._chamber_keys: # First go through chamber products
      try:
        self._add_element(self[key], data.prod_c[key])
      except KeyError:
        self[key].append(0)
    for key in self._exit_keys: # Then go through exit products
      try:
        self._add_element(self[key], data.prod_e[key])
      except KeyError:
        self._add_element(self[key], 0)

if __name__ == "__main__":
  # Example: Set up a simple Water + Peroxide reaction
  water = Fuel("H2O(L)")
  hPeroxide = Oxidizer("H2O2(L)")
  
  chamber_pressure = 1000
  exit_pressure = 14.7
  
  problem = Rocket_Problem(pressure=1000, o_f=5)
  
  for ae_at in [chamber_pressure/exit_pressure, 80]:
    problem.set_ae_at(ae_at)
    data = problem.run_cea(water, hPeroxide)
    print(data)