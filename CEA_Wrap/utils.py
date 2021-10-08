import importlib.resources, os.path

class Output(dict): # This is just a dictionary that you can also use dot notation to access
  def __init__(self): # Explicitly must receive no arguments, because I don't want to deal with constructor properties
    super().__init__()
    
  def __getattr__(self, name):
    return self[name]
  def __setattr__(self, name, value):
    if name.startswith("_"):
      super().__setattr__(name, value)
    else:
      self[name] = value

def _get_asset(file):
  # The reason the manager is used is because our package may be zipped and the manager extracts it
  #   However, this package is not zip-safe so we just return the location
  with importlib.resources.path(__package__+".assets", file) as manager:
    return str(manager)

def open_thermo_lib():
  """
    Opens the attached thermo library input file using the user's default .inp file viewer (should prompt if none)
  """
  print("Opening user manuals using default .inp file viewer")
  os.system(_get_asset("thermo_spg.inp"))

def open_pdfs():
  """
    Opens the attached NASA pdfs using the user's default pdf viewer
  """
  print("Opening user manuals using default pdf viewer")
  os.system('"'+_get_asset("CEA_Mathematical_Analysis.pdf")+'"')
  os.system('"'+_get_asset("CEA_Users_Manual_and_Program_Description.pdf")+'"')
  
def print_assets_directory():
  """
    Just prints the directory where resources are
  """
  var = os.path.dirname(_get_asset("FCEA2.exe"))
  print(var)
  return var