## `function_name(param1: type, param2: type = default, *args, **kwargs) -> return_type`

Short, concise description of what the function does in one line.

### Parameters

- **`param1`** (`type`): Detailed description of the first parameter. Explain its purpose, any constraints, or important characteristics.
  - Example: What kind of input is expected?
  - Any validation or preprocessing requirements?

- **`param2`** (`type`, *optional*, default=`default_value`): Description of the optional parameter.
  - Explain what happens if no value is provided
  - Any specific behaviors or edge cases to be aware of

- **`*args`** (`type`, *optional*): Description of variable positional arguments, if applicable.
  - What kind of additional arguments can be passed?
  - How are they used within the function?

- **`**kwargs`** (`type`, *optional*): Description of variable keyword arguments, if applicable.
  - What additional keyword arguments are supported?
  - How do they modify the function's behavior?

### Returns

- `return_type`: Detailed description of what the function returns.
  - What does the returned value represent?
  - Any special conditions or potential return values?

### Examples

```python
# Basic usage example
result = function_name(arg1, arg2)

# Example with optional parameters
result = function_name(arg1, param2=value)

# Example showing different use cases or edge cases
result = function_name(special_arg, *additional_args)
```

### Notes

- Any additional important information
- Performance considerations
- Potential side effects
- Compatibility notes

# TODO FOR VERSION 2.0

Change how CEA is called to call with location of thermo.lib, location of trans.lib, input is the .inp file, output is the .out/.plt file
  Upon instantiation, we should move our copy of the thermo.lib and trans.lib to a temporary folder (this is for multiprocessing/multithreading)
    Perhaps we only do this if multiprocessing/multithreading is detected?
Change how we actually parse CEA to get the vast majority of things from .plt file instead
Add hecking types to everything! Make CEA output a dataclass. Geez.

Add docstrings to all major functions

## Environment Variables

CEA_Wrap offers a variety of environment variables to customize the locations and behavior when accessing necessary files. In the below options, a "Truthy" value is something like "true", "yes", or "1". A "Falsy" value is something like "false", "no", or "0". This determination only looks at the first character and is case insensitive.

- CEA_USE_SITE_PACKAGES: If set to some Truthy value, CEA_Wrap will not copy assets files to a local directory, instead using the site-packages where CEA_Wrap is installed. This may desired if using CEA_Wrap in a portable installation. This is not recommended if you make modifications to the thermo_spg.inp file, thermo.lib, or trans.lib, because these files will be overwritten if you update CEA_Wrap. DEFAULT: False
- CEA_ASSETS_DIR: If set, CEA_Wrap will use this directory as the "local" directory for finding "assets" such as the CEA executable, thermo.lib, and documentation PDFs. If this is set, the "appdirs" library dependency is not required. This is ignored if CEA_USE_SITE_PACKAGES is set. DEFAULT: None (uses appdirs.user_data_directory)
- CEA_EXE_LOCATION: If set, overrides the path to the CEA executable used for running CEA calculations
- CEA_THERMO_LIB: If set, overrides the path to the thermo.lib file used by the CEA executable
- CEA_TRANS_LIB: If set, overrides the path to the trans.lib file used by the CEA executable
- CEA_THERMO_INP: If set, overrides the path to the thermo_spg.inp file used by the ThermoInterface in CEA_Wrap
- CEA_COPY_THERMO_TO_TEMP: If set to some Truthy value, when a CEA object is instantiated (such as when you import CEA_Wrap), the thermo.lib and trans.lib from assets are copied to a temporary directory (by default). This enables multiprocessing because CEA holds a lock on the thermo.lib and trans.lib files while it is running, and will crash the program if two instances of CEA attempt to access the thermo.lib files at the same time. DEFAULT: True
- CEA_USE_LEGACY: Uses a modified CEA interface that uses legacy logic for interacting with the original FCEA2 executable

# CEA_Wrap

A Python-Based wrapper for the NASA CEA Thermochemical Code

# Installation Instructions

We are now on [PyPi!](https://pypi.org/project/CEA-Wrap/)
In a command prompt type ```pip install --upgrade CEA_Wrap``` to upgrade/install CEA_Wrap

You can import it as any other python module with ```import CEA_Wrap```. Whenever you import the file, it will put the required thermo.lib and trans.lib files into your current directory.

## Installation on Mac/Linux

On mac, the installation script will attempt to compile the fortran executable on your system. Should it fail to do so, you will have to compile it manually.
On Linux, you will need to compile the fortran executable yourself.
You can see a discussion from a successful user here: [Issue #1](https://github.com/civilwargeeky/CEA_Wrap/issues/1#issuecomment-1033918162)

## Making A Portable Installation

As of version 1.7.4 (commit d4331d7), CEA Wrap can be used in a portable manner (without accessing user's program files or app data folders for storing assets).
To download the package to be used as a portable installation, simply clone/download this repository.

To use this package in a portable manner, you must set the environment variable "CEA_ASSETS_DIR" to the directory where "assets" can be found **before** you import CEA_Wrap.
Example structure:

```
Your Code
│   my_code.py 
└───CEA_Wrap
│   │   __init__.py
│   │   __main__.py
│   │   CEA.py
│   │   ...
│   └───assets
│       │   FCEA2.exe
│       │   thermo.lib
│       │   ...
```

Then in my_code.py you would do something like

```python
import os
os.environ["CEA_ASSETS_DIR"] = os.path.join("CEA_Wrap", "assets")
from CEA_Wrap import Oxidizer, Fuel, RocketProblem
...
```

# Examples

Examples on basic use can be found in the "examples" directory. Feel free to download them and try them out!

Or you can run a short demo by doing "python -m CEA_Wrap" on the command line!

# Documentation

### Material(name: str, temp: float = 298.15, wt_percent: float = None, mols: float = None, chemical_composition: str = None, hf: float = None) -> Material

Creates a new material object with specified parameters. All parameters are also members of the class that have both getters and setters

#### Parameters

- **`name`** (`str`): The CEA material name with correct spelling. For example, aluminum is "AL(cr)" and methane is "CH4".
  - If you specify `chemical_composition`, the `name` can be any single word.
  - If you try to specify a chemical which is not in the thermo_spg.inp file, a ValueError will be raised. To prevent this check, set `Material.check_against_thermo_inp = False`.

- **`temp`** (`float`, _optional_, default=`298.15`): The specified reactant initial temperature in Kelvin.
  - If you try to specify a temperature which is not supported for this material in the thermo_spg.inp file, a ValueError will be raised. To prevent this check, set `Material.check_against_thermo_inp = False`.

- **`wt_percent`** (`float`, _optional_, default=`None`): A weight-based percentage for the element.
  - Weight percentages do not need to add up to 100 and are calculated on the ratio with other Fuels/Oxidizers.

- **`mols`** (`float`, _optional_, default=`None`): A mol-based percentage for the element. Can be used as in `Oxidizer("O2", mols=1)` and `Oxidizer("N2", mols=3.76)` for air. (Note: CEA has "Air" as a separate reactant)
  - **NOTE:** `wt_percent` and `mols` cannot be specified together. If neither is defined, the `Material` gets a `wt_percent` of 1.

- **`chemical_composition`** (`str`, _optional_, default=`None`): Chemical composition such as "LI 1 B 1 H 4" for LiBH4.
  - If defined, it will not use CEA default values. **NOTE: Not rigorously tested**

- **`hf`** (`float`, _optional_, default=`None`): Enthalpy of formation in kJ/mol. Must be specified if `chemical_composition` is specified.
  - **NOTE: Not rigorously tested**

#### Returns

- `Material`: A new material object with the specified parameters.

#### Examples

```python
# Basic usage example
material = Material("CH4")

# Example with optional parameters
material = Material(name="AL(cr)", temp=300, wt_percent=50)

# Example showing different use cases or edge cases
material = Material(chemical_composition="LI 1 B 1 H 4", hf=-20)
```

#### Notes

- All parameters from the constructor are also members.
- `.ref` - If `Material.check_against_thermo_inp` is True, this will be a `ThermoMaterial` representing the material.

### Available Methods

#### set_wt_percent(wt_percent: float) -> None

Sets the wt_percent for the Material. Sets mols to None.

```python
material.set_wt_percent(50)
```

#### set_mols(mols: float) -> None

Sets the mols for the Material. Sets wt_percent to None.

```python
material.set_mols(1)
```

#### set_temp(temp: float) -> None

Sets the temp for the Material, in Kelvin. If you try to specify a temperature which is not supported for this material, a ValueError will be raised.

```python
material.set_temp(300)
```

#### is_mols() -> bool

Checks if the Material is set in mol ratio.

```python
is_mol = material.is_mols()
```

#### is_wt_percent() -> bool

Checks if the Material is set in weight percentage.

```python
is_wt_percent = material.is_wt_percent()
```

## `Problem(*, **kwargs) -> Problem`

Creates a new Problem object with specified parameters.

### Parameters

- **`**kwargs`** (`dict`, *optional*): Keyword arguments for initializing the problem.
  - All parameters must be specified by keyword, e.g., `problem = RocketProblem(pressure=500, massf=True)`.
  - Example: What kind of input is expected?
    - Any validation or preprocessing requirements?

- **`pressure`** (`float`, *optional*, default=1000): Initial problem pressure.
  - This parameter sets the initial pressure for the problem.

- **`materials`** (`list`, *optional*, default=None): List of `Material` objects, order doesn't matter, of Oxidizer and Fuel objects e.g. materials=[material1, material2, ...].
  - Materials can also be specified when you run a problem like `problem.run_cea(material1, material2, ...)`.
  - Materials MUST all have wt_percent specified or all have mols specified; cannot have mixtures.

- **`massf`** (`bool`, *optional*, default=False): CEA usually outputs product species in mole fractions. If massf is True, mass fractions are specified.
  - This parameter determines whether the output should be in mole fractions or mass fractions.

- **`filename`** (`str`, *optional*, default="my_output"): The filename we save our .inp files to and CEA saves our .out and .plt to.
  - DO NOT INCLUDE ".inp" IN YOUR FILENAMES.

- **`pressure_units`** (`str`, *optional*, default="psi"): The units that your input pressure is in. Possible values are "bar", "atm", "psi", or "mmh".
  - This parameter sets the units for the input pressure.

- **`inserts`** (`list or str`, *optional*, default=None): A list of CEA names for species which should be forced into the product considerations. Should be specified as either a space-separated string of names, a list of string names, or a list of `Material` objects.
  - Note: If you try to specify a chemical which is not in the thermo_spg.inp file, a ValueError will be raised. To prevent this check, set `Problem.check_against_thermo_inp = False`.
  - Tip: If you are doing calculations with Aluminum, I recommend using inserts=["AL2O3(L)", "AL2O3(a)"].

- **`omits`** (`list or str`, *optional*, default=None): A list of CEA names for species which should be specifically ignored in the product considerations. Specified similar to inserts.
  - This parameter allows you to exclude certain species from the product considerations.

#### Specifying reactant ratios:

   Key  | CEA Key | Description
  ------|---------|---------------------------------------------------------------------------------------------------------
   p_f  |  %f     | Percent fuel by weight
   f_o  |  f/o    | Fuel-to-oxidant weight ratio
   o_f  |  o/f    | Oxidant-to-fuel weight ratio
   phi  |  phi    | Equivalence ratios in terms of fuel-to-oxidant weight ratios (eq. (9.19) in Gordon and McBride, 1994)
   r_eq |  r      | Chemical equivalence ratios in terms of valences (eq. (9.18) in Gordon and McBride, 1994)


### Returns

- `Problem`: An instance of the Problem class initialized with the provided parameters.

### Examples

```python
# Basic usage example
problem = Problem(pressure=1000, materials=[material1, material2], massf=True)

# Example with optional parameters
problem = Problem(filename="custom_output", pressure_units="bar")
```

### Notes

- All parameters must be specified by keyword.
- Performance considerations: Initializing a problem object may involve significant computation depending on the complexity of the input parameters.

## `Problem.run([*materials]) -> Output`

Runs the CEA problem and returns an "Output" object similar to a dictionary.

### Parameters

- **`*materials`** (`list`, *optional*): List of `Material` objects to be used in this run.
  - If materials are not specified as an initial parameter or with `.set_materials()`, you can list them here.

### Returns

- `Output`: An "Output" object containing the results of the CEA problem, similar to a dictionary.

### Examples

```python
# Basic usage example
data = problem.run()

# Example with materials
data = problem.run(material1, material2)
```

## `Problem.set_absolute_o_f() -> None`

Calculates and sets the correct o/f ratio based on the existing material list.

### Parameters

- **None**

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_absolute_o_f()
```

## `Problem.set_pressure(pressure) -> None`

Sets the pressure for the problem.

### Parameters

- **`pressure`** (`float`): The new pressure value to set.

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_pressure(100)
```

## `Problem.set_materials([material1, material2, ...]) -> None`

Sets the materials list for the problem.

### Parameters

- **`[material1, material2, ...]`** (`list`): List of `Material` objects to set as the materials list.

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_materials([material1, material2])
```

## `Problem.set_massf(massf) -> None`

Sets the mass fraction flag for the problem.

### Parameters

- **`massf`** (`bool`): The new value for the mass fraction flag (True or False).

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_massf(True)
```

## `Problem.set_inserts(inserts) -> None`

Sets the inserts for the problem.

### Parameters

- **`inserts`** (`str or list`, *optional*): A space-separated string of names, a list of string names, or a list of `Material` objects to set as inserts.
  - If you try to specify a chemical which is not in the thermo_spg.inp file, a ValueError will be raised. To prevent this check, set `Problem.check_against_thermo_inp = False`.

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_inserts(["AL2O3(L)", "AL2O3(a)"])
```

## `Problem.set_omits(omits) -> None`

Sets the omits for the problem.

### Parameters

- **`omits`** (`str or list`, *optional*): A space-separated string of names, a list of string names, or a list of `Material` objects to set as omits.
  - Specified similar to inserts.

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_omits(["CO2", "H2O"])
```

## `Problem.set_filename(filename) -> None`

Sets the filename for the problem.

### Parameters

- **`filename`** (`str`): The new filename to set.
  - DO NOT INCLUDE ".inp" IN YOUR FILENAMES.

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_filename("custom_output")
```

## `Problem.set_pressure_units(units) -> None`

Sets the pressure units for the problem.

### Parameters

- **`units`** (`str`): The new pressure units to set. Possible values are "bar", "atm", "psi", or "mmh".

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_pressure_units("bar")
```

## `Problem.set_p_f(pf) -> None`

Sets the percent fuel by weight for the problem.

### Parameters

- **`pf`** (`float`): The new value for percent fuel by weight.

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_p_f(0.5)
```

## `Problem.set_f_o(f_o) -> None`

Sets the fuel-to-oxidant weight ratio for the problem.

### Parameters

- **`f_o`** (`float`): The new value for the fuel-to-oxidant weight ratio.

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_f_o(2.0)
```

## `Problem.set_o_f(o_f) -> None`

Sets the oxidizer-to-fuel weight ratio for the problem.

### Parameters

- **`o_f`** (`float`): The new value for the oxidizer-to-fuel weight ratio.

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_o_f(0.5)
```

## `Problem.set_phi(phi) -> None`

Sets the equivalence ratio for the problem.

### Parameters

- **`phi`** (`float`): The new value for the equivalence ratio.

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_phi(1.2)
```

## `Problem.set_r_eq(r_eq) -> None`

Sets the chemical equivalence ratio for the problem.

### Parameters

- **`r_eq`** (`float`): The new value for the chemical equivalence ratio.

### Returns

- `None`: This method does not return any value; it modifies the object in place.

### Examples

```python
# Basic usage example
problem.set_r_eq(1.5)
```

## TP (Specified Temperature and Pressure) Problem Constructor Additional Parameters:

Very similar to an HP problem, but temperature is specified per-problem and material temperatures are ignored.

### Parameters

- **`temperature`** (`float`, *optional*, default=298): The problem temperature.
  - Can be specified later with `.set_temperature`.
- **`temperature_units`** (`str`, *optional*, default='k'): The units for the temperature. Options are 'k', 'c', 'r', 'f'.
  - Can be specified later with `.set_temperature_units`.

### Returns

- `TPProblem`: An instance of the TPProblem class initialized with the provided parameters.

## Detonation Problem Constructor Additional Parameters:

**WARNING**: As far as I am aware, CEA is incapable of performing detonation calculations with condensed phase (solid) reactants. It will only handle gaseous reactants.

### Returns

- `DetonationProblem`: An instance of the DetonationProblem class initialized with the provided parameters.

## Rocket Problem Constructor Additional Parameters:

For `RocketProblem(*, **kwargs)`

### Parameters

- **`sup`** (`float`, *optional*, default=1): Supersonic exit area/throat area ratio.
  - Can be specified later with `.set_sup`.
- **`sub`** (`float`, *optional*): Subsonic exit/throat area ratio.
  - Can be specified later with `.set_sub`.
  - **NOTE:** `sup` and `sub` cannot be specified at the same time.
- **`ae_at`** (`float`, *optional*, alias for `sup`): Alias for supersonic exit area/throat area ratio.
  - Can be specified later with `.set_ae_at`.
- **`pip`** (`float`, *optional*): Pressure ratio of chamber pressure/exit pressure.
  - Can be specified later with `.set_pip`.
- **`fac_ma`** (`float`, *optional*): Finite Area Combustor, with mass flow (mdot) / combustor chamber area. Units of (kg/s)/m^2.
  - Can be specified later with `.set_fac_ma`.
  - Cannot be specified at the same time as `fac_ac`.
- **`fac_ac`** (`float`, *optional*): Finite Area Combustor, with ac/at: Ratio of combustor area to throat area.
  - Can be specified later with `.set_fac_ac`.
  - Cannot be specified at the same time as `fac_ma`.
- **`analysis_type`** (`str`, *optional*, default="equilibrium"): Whether to use equilibrium reactions or frozen. For using frozen specify "frozen" or "frozen nfz=1" for frozen at the chamber or "frozen nfz=2" for frozen at the throat.
- **`nfz`** (`int`, *optional*): If `analysis_type` is "frozen", this will set the frozen location to the given point. 1 for chamber, 2 for throat.
- **`custom_nfz`** (`float`, *optional*): If `analysis_type` is "frozen", this is the position within the nozzle that composition will be frozen at. Uses the same unit as ae/at or pip.

### Returns

- `RocketProblem`: An instance of the RocketProblem class initialized with the provided parameters.

## Using ThermoInterface

Provided with this library is an interface to the thermo_spg.inp file provided with the library.

You can access materials as if `ThermoInterface` is a dictionary using materials' CEA Names. E.G. `ThermoInterface["CH4"]`. This object support checking for inclusion and iterations, such as `"CH4" in ThermoInterface` and `for name in ThermoInterface`. It supports dictionary methods such as `.keys()`, `.values()` and `.items()`.

The value returned by `ThermoInterface` accesses is a `ThermoMaterial object`, which is an `Output` object (dictionary that can be accessed with .dot notation) with the following keys:

- `name` - Name of the material
- `reference` - Reference, if given. Otherwise ""
- `elements` - Dictionary of element: numerical value specified
- `condensed` - True if condensed phase, False otherwise
- `mw` - Molecular weight in g/mol
- `hf` - float heat of formation at 298.15 in kJ/mol (or assigned enthalpy if 0 temp range)
- `temp_ranges` - List of two-tuples of [range start, range end] (K)
- `reactant_only` - True if material only shows up in reactants, False otherwise

### Available Methods

- `ThermoInterface.get_close_matches(name, [n])` - Gets close matches to a given material name. For example "Al(cr)" returns 'AL(cr)', 'ALN(cr)', 'Ag(cr)', 'W(cr)'. n influences the number of results returned, and is the maximum number of results returned.
- `ThermoMaterial.defined_at(temp)` - Returns True if the material is specified at the given temperature, False otherwise. Materials specified at one temperature are actually allowed at that temperature +- 10K.

## Utilities

  ```open_thermo_lib()```
  Opens the default thermo library input file using the user's default .inp file viewer (should prompt if none)
  
  ```open_pdfs()```
  Opens the attached NASA pdfs using the user's default pdf viewer
  
  ```print_assets_directory()```
  Prints to console the current location of the directory where CEA_Wrap assets are located. Also returns this value
  
  ```print_simple_thermo_lib_line(name, comment, atoms, isCond, molWt, hf, temperature=298.15)```
  Returns and prints a string which can be inserted to represent a molecule in the Thermo Lib. Any entries which are longer than the allotted space results in an error
  
- `name`: 24 chars max, Species name, such as CO or Fe2O3. These are assumed to be in a gas phase unless appended with (a) for agglomerate, (cr) for crystalline (solid), or (L) for liquid phase.
  
- `comment`: 56 char max, Human-readable name and additional information
  
- `atoms`: 5 atom max, A dictionary of "atomic symbol": number of atoms in molecule.
  
    Example: Water, "H2O", would be {"H": 2, "O": 1}

    Complex Example: Lead Acetate, "Pb(C2H3O2)2", would be {"Pb": 1, "C": 4, "H": 6, "O": 4}. Note that the C, H, and O are doubled because they are the *sum* of atoms in the molecule

    Note: "E" is a special value which represents an electron for ionic compounds. You can have negative amounts of "E" to indicate positively charged atoms

- `isCond`: "Is Condensed?" True if the material is a solid or liquid, False if it is a gas.
  
- `molWt`: Molecular weight of molecule in g/mol or kg/kmol
  
- `hf`: Assigned Enthalpy of material at temperature specified. J/mol. If temperature=298.15K, this is the Heat of Formation.
  
- `temperature`: Temperature that the molecule is specified at, Kelvin
  
  ```reload_thermo_lib()```
  Moves the thermo and transport libs from your assets directory (if they have changed) and reloads ThermoInterface from the thermo_spg.inp in your assets directory
  Use this if you have recompiled your thermo lib but do not want to restart python.
  
  Note: This is actually defined in CEA.py, not utils.py
  
### DataCollector

  A DataCollector object will conveniently compile output from Problem outputs into list format. Example usage can be found in the "looping.py" example
  
  ```DataCollector(self, *args, keys=[], chamber_keys=[], throat_keys=[], exit_keys=[])```
  
- `keys`: Also accepts list of arguments, these are keys such as 'cond' or 't_cp' or 'c_p'
- `chamber_keys`: List of chamber species mol/mass fractions to be included in the output object. Ex: "H2O" or "CO2". The key in the output will be the molecule name with "c_" prepended
- `throat_keys`: List of nozzle throat species mol/mass fractions to be included. The key in the output will be the molecule name with "t_" prepended
- `exit_keys`: List of nozzle exit species mol/mass fractions to be included. The key in the output will be the molecule name with nothing prepended.
  
  *Methods*
- `add_data(data)` - Data should be the output from Problem.run(). Appends the output of the keys specified in the initializer to the object
- `to_csv(filename, filename: str, keys:list=None, formatString="f"` - Writes the data to csv at **filename**, with the keys in **keys** (or all keys in this object if None)
     using **formatString** to format the csv entries

## Adding materials/recompiling the thermo lib (advanced usage)
  
  You can add new molecules by modifying your thermo_spg.inp file and then recompiling the thermo lib.

  1. Find the location of your thermo_spg.inp file by using print_assets_directory()
  1. Open thermo_spg.inp in a text editor.
  1. Scroll to the bottom
  1. You can use print_simple_thermo_lib_line() to create a simple entry in the format of other reactants
  1. Add your entry below the other entries but above "END REACTANTS". Ensure your entry has similar spacing to existing entries
  1. To recompile your thermo lib, either use the "thermo_lib_recompile" batch file or run "FCEA2" and type "thermo_spg"
  1. If you have a python kernel open, you should now restart it or use the "reload_thermo_lib()" function to reload the thermo lib used.
