# CEA_Wrap

A Python-Based wrapper for the NASA CEA Thermochemical Code

# Installation Instructions

We are now on [PyPi!](https://pypi.org/project/CEA-Wrap/)
In a command prompt type ```pip install --upgrade CEA_Wrap``` to upgrade/install CEA_Wrap

You can import it as any other python module with ```import CEA_Wrap```

## Installation on Mac/Linux

On mac, the installation script will attempt to compile the fortran executable on your system. Should it fail to do so, you will have to compile it manually.
On Linux, you will need to compile the fortran executable yourself.
You can see a discussion from a successful user here: [Issue #1](https://github.com/civilwargeeky/CEA_Wrap/issues/1#issuecomment-1033918162)

## Portable Installation

As of version 2.0, making a portable installation is even easier! This will initialize CEA to set itself up without accessing user's program files or app data folders for storing assets. Simply set the `CEA_USE_SITE_PACKAGES` environment variable before loading the code.

Details can be found in "Making A Portable Installation" at the bottom of the README

# Important: Changes to transport properties in Version 2.1

It was brought to my attention through an email exchange and [Issue #11](https://github.com/civilwargeeky/CEA_Wrap/issues/11) that transport properties in CEA_Wrap prior to version 2.1 do not match those from the CEARun website. This is because the "trans.lib" file I was originally handed down did not use the coefficients described in the CEA documentation, the trans.lib I had was blank... This means that calculations were not using experimentally-derived transport properties, and were using less-accurate fallback methods.

I **strongly urge you to update** to `pip install --upgrade CEA_Wrap` to a version >= 2.1 to get more-accurate transport properties in your calculations.

# New in Version 2

CEA_Wrap has undergone a partial rewrite in Version 2.0 to change the underlying way that many of the mechanisms work and more type hinting for those using IDEs to interact with it. A summary of major changes:

1. **Enabled Multithreaded Support**: By creating a temporary directory containing the thermo.lib/trans.lib for each thread/process, we allow for hugely increased processing speeds when using multiple threads
1. Vastly improved function documentation and type-hinting. For example: outputs from Problem.run() are now dataclasses, so their members are type-hinted. All major functions should have docstrings, so that cross-referencing this README is less-needed
1. Rewrote some of cea.f to accept input from STDIN and write output to STDOUT. This eliminates the need to write/read .inp/.out/.plt files and speeds up processing on slower filesystems
1. Added a variety of environment variables to allow for customization of interaction with CEA. Allows for using multiple thermo.lib files, for example.

# Examples

Examples on basic use can be found in the "examples" directory. Feel free to download them and try them out!

Or you can run a short demo by doing "python -m CEA_Wrap" on the command line!

# Documentation

## Defining Materials

### class `Material(name: str, temp: float = 298.15, wt: float|None = None, mols: float|None = None, chemical_representation: ChemicalRepresentation|None = None, hf: float|None = None) -> Material`

Creates a new material object with specified parameters. All parameters are also members of the class that have both getters and setters. You should specify either `wt` or `mols` but cannot specify both. If neither is defined, `wt` is set to `1`.

##### Parameters

- **`name`** (`str`): The CEA material name with correct spelling. For example, aluminum is "AL(cr)" and methane is "CH4".
  - If you specify `chemical_composition`, the `name` can be any single word.
  - If you try to specify a chemical which is not in the thermo_spg.inp file, a ValueError will be raised. To prevent this check, set `Material.check_against_thermo_inp = False`.

- **`temp`** (`float`, _optional_, default=`298.15`): The specified reactant initial temperature in Kelvin.
  - If you try to specify a temperature which is not supported for this material in the thermo_spg.inp file, a ValueError will be raised. To prevent this check, set `Material.check_against_thermo_inp = False`.

- **`wt`** (`float`, _optional_, default=`None`): A weight-based percentage for the element.
  - Weight percentages do not need to add up to 100 and are calculated with ratios according to other Fuels/Oxidizers in a given problem.

- **`mols`** (`float`, _optional_, default=`None`): A mol-based percentage for the element. Can be used as in `Oxidizer("O2", mols=1)` and `Oxidizer("N2", mols=3.76)` for air. (Note: CEA has "Air" as a separate reactant)

- **`chemical_representation`** (`ChemicalRepresentation`, _optional_, default=`None`): Chemical composition such as "LI 1 B 1 H 4" for LiBH4.
  - If defined, it will not use CEA default values. **NOTE: Not rigorously tested**

##### Examples

```python
# Basic usage example
material = Material("CH4")

# Example with optional parameters
material = Material(name="AL(cr)", temp=300, wt=50)

# Example with specifying chemical_composition
material = Material(chemical_composition="LI 1 B 1 H 4", hf=-20)
```

##### Notes

- All parameters from the constructor are also members.
- `.ref` - If `Material.check_against_thermo_inp` is True, this will be a `ThermoMaterial` representing the material.

#### Available Methods

##### set_wt_percent(wt_percent: float) -> None

Sets the wt_percent for the Material. Sets mols to None.

##### set_mols(mols: float) -> None

Sets the mols for the Material. Sets wt_percent to None.

##### set_temp(temp: float) -> None

Sets the temp for the Material, in Kelvin. If you try to specify a temperature which is not supported for this material, a ValueError will be raised.

##### is_mols() -> bool

Checks if the Material is set to represent a number of mols.

##### is_wt() -> bool

Checks if the Material is set to represent weight percentage / mass.

### dataclass `ChemicalRepresentation(chemical_composition: str, hf: float, hf_unit: str = "kj", is_internal_energy:bool = False)`
A ChemicalRepresentation represents a custom material added in the .inp file, as opposed to using an existing material in thermo.lib

##### Parameters

- **`chemical_representation`** (`str`): chemical composition in space-separated pairs of (element, n_atoms), such as "LI 1 B 1 H 4" for LiBH4.
  
- **`hf`** (`float`): Enthalpy of formation, in `hf_unit`/mol.

- **`hf_unit`** (`float`, _optional_, default=`kj`): Unit for enthalpy of formation.
  - `hf_unit` can be any of "c" (calorie), "kc" (kilo-calorie), "j" (joule), or "kj" (kilo-joule)

- **`is_internal_energy`** (`bool`, _optional_, default=`False`): If True, sets an internal energy (`u`), rather than enthalpy (`h`).

### `Problem(*, **kwargs) -> Problem`

Creates a new Problem object with specified parameters.

#### Parameters

- **`**kwargs`** : Keyword arguments for initializing the problem.
  - All parameters must be specified by keyword, e.g., `problem = RocketProblem(pressure=500, massf=True)`.

- **`pressure`** (`float`, *optional*, default=1000): Initial problem pressure.
  - This parameter sets the initial pressure for the problem.

- **`materials`** (`list`, *optional*, default=None): List of `Material` objects, order doesn't matter, of Oxidizer and Fuel objects e.g. materials=[material1, material2, ...].
  - Materials can also be specified when you run a problem like `problem.run_cea(material1, material2, ...)`.
  - Materials MUST all have wt specified or all have mols specified; cannot have mixtures.
  - You can modify a Material object, for instance changing it's `wt`, and re-run the problem. You do not need to set the materials again after changing weights.

- **`massf`** (`bool`, *optional*, default=False): CEA usually outputs product species in mole fractions. If massf is True, mass fractions are specified.
  - This parameter determines whether the output should be in mole fractions or mass fractions.

- **`pressure_units`** (`str`, *optional*, default="psi"): The units that your input pressure is in. Possible values are "bar", "atm", "psi", or "mmh".
  - This parameter sets the units for the input pressure.

- **`inserts`** (`list or str`, *optional*, default=None): A list of CEA names for species which should be forced into the product considerations. Should be specified as either a space-separated string of names, a list of string names, or a list of `Material` objects.
  - Note: If you try to specify a chemical which is not in the thermo_spg.inp file, a ValueError will be raised. To prevent this check, set `Problem.check_against_thermo_inp = False`.
  - Tip: If you are doing calculations with Aluminum, I recommend using inserts=["AL2O3(L)", "AL2O3(a)"]. This may help errors with convergence at lower temperatures/pressures. You can also do similar things for other condensed phase species.

- **`omits`** (`list or str`, *optional*, default=None): A list of CEA names for species which should be specifically ignored in the product considerations. Specified similar to inserts.
  - This parameter allows you to exclude certain species from the product considerations.

##### Specifying reactant ratios:

   Key  | CEA Key | Description
  ------|---------|---------------------------------------------------------------------------------------------------------
   p_f  |  %f     | Percent fuel by weight
   f_o  |  f/o    | Fuel-to-oxidant weight ratio
   o_f  |  o/f    | Oxidant-to-fuel weight ratio
   phi  |  phi    | Equivalence ratios in terms of fuel-to-oxidant weight ratios (eq. (9.19) in Gordon and McBride, 1994)
   r_eq |  r      | Chemical equivalence ratios in terms of valences (eq. (9.18) in Gordon and McBride, 1994)


#### Returns

- `Problem`: An instance of the Problem class initialized with the provided parameters.
  - Note: The base `Problem` class cannot be run by itself, you must run a sub-class such as HPProblem or RocketProblem

#### Examples

```python
# Basic usage example
problem = Problem(pressure=1000, materials=[material1, material2], massf=True)

# Example with optional parameters
problem = Problem(inserts=["AL2O3(L)", "AL2O3(a)"], pressure_units="bar")
```

#### Notes

- All parameters must be specified by keyword.

### `Problem.run([*materials]) -> Output`

Runs the CEA problem and returns an "Output" object similar to a dictionary.

#### Parameters

- **`*materials`** (`list`, *optional*): List of `Material` objects to be used in this run.
  - If materials are not specified as an initial parameter or with `.set_materials()`, you can list them here.

#### Returns

- `Output`: An "Output" object containing the results of the CEA problem, similar to a dictionary. Elements of the Output may be accessed with ["braces"] or .dot_notation
  - Example: `problem.run()["T_c"]` or `problem.run().T_c`
  - As of Version 2.0, all Problem subclasses return a problem-specific dataclass type, exposing all the parameters available for autocomplete in editors like VS Code or PyCharm.

#### Examples

```python
# Basic usage example
data = problem.run()

# Example with materials
data = problem.run(material1, material2)
```

### `Problem.set_absolute_o_f() -> None`

Calculates and sets the correct o/f ratio based on the existing material list.

#### Parameters

- **None**

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem = Problem(materials=[Fuel("A", wt=30), Oxidizer("B", wt=60)], o_f=1)
problem.set_absolute_o_f() # Sets o_f to 60/30
```

### `Problem.set_pressure(pressure) -> None`

Sets the pressure for the problem.

#### Parameters

- **`pressure`** (`float`): The new pressure value to set.

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem.set_pressure(100)
```

### `Problem.set_materials([material1, material2, ...]) -> None`

Sets the materials list for the problem.

#### Parameters

- **`[material1, material2, ...]`** (`list`): List of `Material` objects to set as the materials list.

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem.set_materials([material1, material2])
```

### `Problem.set_massf(massf) -> None`

Sets the mass fraction flag for the problem. If massf is True, CEA returns mass-fractions for exhaust product abundances, otherwise it returns mol-fractions (the default)

#### Parameters

- **`massf`** (`bool`): The new value for the mass fraction flag (True or False).

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem.set_massf(True)
```

### `Problem.set_inserts(inserts) -> None`

Sets the inserts for the problem.

#### Parameters

- **`inserts`** (`str or list`, *optional*): A space-separated string of names, a list of string names, or a list of `Material` objects to set as inserts.
  - If you try to specify a chemical which is not in the thermo_spg.inp file, a ValueError will be raised. To prevent this check, set `Problem.check_against_thermo_inp = False`.

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem.set_inserts(["AL2O3(L)", "AL2O3(a)"])
```

### `Problem.set_omits(omits) -> None`

Sets the omits for the problem.

#### Parameters

- **`omits`** (`str or list`, *optional*): A space-separated string of names, a list of string names, or a list of `Material` objects to set as omits.
  - Specified similar to inserts.

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem.set_omits(["CO2", "H2O"])
```

### `Problem.set_pressure_units(units) -> None`

Sets the pressure units for the problem.

#### Parameters

- **`units`** (`str`): The new pressure units to set. Possible values are "bar", "atm" (atmospheres), "psi" (pounds per square inch), or "mmh" (millimeters mercury).

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem.set_pressure_units("bar")
```

### `Problem.set_p_f(pf) -> None`

Sets the percent fuel by weight for the problem.
**WARNING**: There seems to be a bug in CEA itself that this option is simply interpreted as an alias for O/F ratio

#### Parameters

- **`pf`** (`float`): The new value for percent fuel by weight.

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem.set_p_f(0.5)
```

### `Problem.set_f_o(f_o) -> None`

Sets the fuel-to-oxidant weight ratio for the problem.

#### Parameters

- **`f_o`** (`float`): The new value for the fuel-to-oxidant weight ratio.

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem.set_f_o(2.0)
```

### `Problem.set_o_f(o_f) -> None`

Sets the oxidizer-to-fuel weight ratio for the problem.

#### Parameters

- **`o_f`** (`float`): The new value for the oxidizer-to-fuel weight ratio.

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem.set_o_f(0.5)
```

### `Problem.set_phi(phi) -> None`

Sets the equivalence ratio for the problem, using the engineering definition of equivalence ratio.

#### Parameters

- **`phi`** (`float`): The new value for the equivalence ratio.

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem.set_phi(1.2) # Fuel-rich
```

### `Problem.set_r_eq(r_eq) -> None`

Sets the 'chemical equivalence ratio' for the problem.

#### Parameters

- **`r_eq`** (`float`): The new value for the chemical equivalence ratio.

#### Returns

- `None`: This method does not return any value; it modifies the object in place.

#### Examples

```python
# Basic usage example
problem.set_r_eq(1.5)
```

### HP (Specified Enthalpy and Pressure) Problem Constructor Additional Parameters:

Usage: `HPProblem(**kwargs)`.

This problem type has no additional parameters beyond those in `Problem`. In an HP problem, the enthalpy of the system is calculated based on the temperatures of the individual reactant materials specified.

### TP (Specified Temperature and Pressure) Problem Constructor Additional Parameters:

Very similar to an HP problem, but temperature is specified per-problem and thus material temperatures are ignored.

#### Parameters

- **`temperature`** (`float`, *optional*, default=298): The problem temperature.
  - Can be specified later with `.set_temperature`.
- **`temperature_units`** (`str`, *optional*, default='k'): The units for the temperature. Options are 'k' (Kelvin), 'c' (Celsius/Centigrade), 'r' (Rankine), 'f' (Fahrenheit).
  - Can be specified later with `.set_temperature_units`.

#### Returns

- `TPProblem`: An instance of the TPProblem class initialized with the provided parameters.

### Detonation Problem Constructor Additional Parameters:

Usage: `DetonationProblem(**kwargs)`. This problem type has no additional parameters beyond those in `Problem`

**NOTICE**: As far as I am aware, CEA is incapable of performing detonation calculations with condensed phase (solid) reactants. It will only handle gaseous reactants.

#### Returns

- `DetonationProblem`: An instance of the DetonationProblem class initialized with the provided parameters.

### Rocket Problem Constructor Additional Parameters:

For `RocketProblem(*, **kwargs)`. The `RocketProblem` has a variety of options related to how the chamber and nozzle of the rocket motor are treated.

#### Parameters

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
- **`custom_nfz`** (`float`, *optional*): If `analysis_type` is "frozen", this is the position within the nozzle that composition will be frozen at. Uses the same unit as ae/at or pip (whichever is set).

#### Returns

- `RocketProblem`: An instance of the RocketProblem class initialized with the provided parameters.

## Available Output Dictionary Keys:
All Problem data objects are "Output" objects, which are similar to dictionaries, but can also be accessed with dot notation.

For example if you had "data = problem.run_cea()", and wanted pressure, you could do either data.p or data["p"]

In addition, all product dictionaries are also "Output" objects so to get H2O fraction, you could use data.prod_c.H2O or data["prod_c"]["H2O"] or data["prod_c"].H2O, etc.

### Detonation:
* `prod_c` - dictionary of chamber products, in mole or mass fractions (as specified by the 'massf' parameter in your problem definition)
* `massf` - True if this problem returns mass fractions, False otherwise
* `p` - pressure, bar
* `t` - temperature, Kelvin
* `h` - enthalpy, kJ/kg
* `rho` - density, kg/m^3
* `son` - sonic velocity, m/s
* `visc` - burned gas viscosity, Pascal-Seconds
* `mw` - molecular weight of products, kg/kmol
* `cp` - constant-pressure specific heat capacity, kJ/(kg*K)
* `gammas` - isentropic ratio of specific heats
* `gamma` - "real" ratio of specific heats (multiplied by -(dLV/dLP)t)
* `vel` - detonation velocity, m/s
* `mach` - detonation mach number
* `p_p1` - P/P1, ratio of detonation pressure to initial pressure
* `t_t1` - T/T1, ratio of detonation temperature to initial pressure
* `m_m1` - M/M1, ratio of detonation molecular weight to initial molecular weight
* `rho_rho1` - RHO/RHO1, ratio of detonation density to initial density
* `dLV_dLP_t` - (dLV/dLP)t, used to convert isentropic gamma to real gamma
* `dLV_dLT_p` - (dLV/dLT)p
* `cond` - burned gas thermal conductivity, W/(m K)
* `pran` - burned gas Prandtl number. Ratio of mass diffusivity to thermal diffusivity
* `phi` - weight-based equivalence ratio of oxidizer/fuel in original problem
### HP (Specified Enthalpy and Pressure) and TP (Specified Temperature and Pressure):
* `prod_c` - dictionary of chamber products, in mole or mass fractions  (as specified by the 'massf' parameter in your problem definition)
* `p` - pressure, bar
* `t` - temperature, Kelvin
* `h` - enthalpy, kJ/kg
* `rho` - density, kg/m^3
* `son` - sonic velocity, m/s
* `visc` - burned gas viscosity, Pascal-Seconds
* `mw` - molecular weight of products, kg/kmol
* `cp` - constant-pressure specific heat capacity, kJ/(kg*K)
* `gammas` - isentropic ratio of specific heats
* `gamma` - "real" ratio of specific heats (multiplied by -(dLV/dLP)t)
* `dLV_dLP_t` - (dLV/dLP)t, used to convert isentropic gamma to real gamma
* `dLV_dLT_p` - (dLV/dLT)p
* `cond` - burned gas thermal conductivity, W/(m K)
* `pran` - burned gas Prandtl number. Ratio of mass diffusivity to thermal diffusivity
* `phi` - weight-based equivalence ratio of oxidizer/fuel in original problem
### Rocket:
* **NOTE : Properties are by default at exit. Chamber parameters are prefixed "c_" and throat properties "t_"**
* **NOTE : Properties not defined for frozen flow are marked with an asterisk (*)**
* **NOTE : Properties not defined _only_ for frozen flow are marked with a dagger (†)**
* **NOTE : All properties defined at the throat are also defined as "f_property" when Finite Area Combustor is enabled (defined fac_ac or fac_ma)**
* `prod_c` - dictionary of chamber products, in mole or mass fractions (as specified by the 'massf' parameter in your problem definition)
* `*prod_t` - dictionary of throat products, in mole or mass fractions (as specified by the 'massf' parameter in your problem definition)
* `*prod_e` - dictionary of exit products, in mole or mass fractions (as specified by the 'massf' parameter in your problem definition)
* `p` - pressure, bar
  * `t_p` - throat
  * `c_p` - chamber
* `t` - temperature, Kelvin
  * `t_t` - throat
  * `c_t` - chamber
* `h` - enthalpy, kJ/kg
  * `t_h` - throat
  * `c_h` - chamber
* `rho` - density, kg/m^3
  * `t_rho` - throat
  * `c_rho` - chamber
* `son` - sonic velocity, m/s
  * `t_son` - throat
  * `c_son` - chamber
* `visc` - burned gas viscosity, Pascal-Seconds
  * `t_visc` - throat
  * `c_visc` - chamber
* `cond` - burned gas thermal conductivity, W/(m K)
  * `t_cond` - throat
  * `c_cond` - chamber
* `†cond_eq` - equilibrium calculated burned gas thermal conductivity, W/(m K) (in frozen flow, most properties represent frozen conditions)
  * `†t_cond_eq` - throat
  * `†c_cond_eq` - chamber
* `pran` - burned gas Prandtl number. Ratio of mass diffusivity to thermal diffusivity
  * `t_pran` - throat
  * `c_pran` - chamber
* `mw` - molecular weight of all products, kg/kmol
  * `t_mw` - throat
  * `c_mw` - chamber
* `m` - molecular weight calculated as the weight of all products divided by the number of gaseous moles (same as mw if no condensed phases as mw=m), kg/kmol
  * `t_m` - throat
  * `c_m` - chamber
* `condensed` - True if there are condensed phase products (measured by m originally = 0), False otherwise
  * `t_condensed` - throat
  * `c_condensed` - chamber
* `cp` - constant-pressure specific heat capacity, kJ/(kg*K)
  * `t_cp` - throat
  * `c_cp` - chamber
* `gammas` - isentropic ratio of specific heats
  * `t_gammas` - throat
  * `c_gammas` - chamber
* `gamma` - "real" ratio of specific heats (multiplied by -(dLV/dLP)t) (same as gammas for frozen flow)
  * `t_gamma` - throat
  * `c_gamma` - chamber
* `isp` - ideal isp (ambient pressure = exit pressure), s
  * `t_isp` - throat
* `ivac` - vacuum isp, s
  * `t_ivac` - throat
* `cf` - ideally expanded thrust coefficient
  * `t_cf` - throat
* `*dLV_dLP_t` - (dLV/dLP)t, multiply gammas by negative this to convert isentropic gamma to real gamma
  * `*t_dLV_dLP_t` - throat
  * `*c_dLV_dLP_t` - chamber
* `*dLV_dLT_p` - (dLV/dLT)p, see the Mathematical Analysis in CEA_Wrap/assets for explanation
  * `*t_dLV_dLT_p` - throat
  * `*c_dLV_dLT_p` - chamber
* `cstar` - characteristic velocity in chamber, m/s
* `*mach` - mach number at exhaust
* `o_f` - oxidizer/fuel weight ratio
* `phi` - weight-based equivalence ratio of oxidizer/fuel
* `ae` - Ratio of area at exit to area at throat.
  * `t_ae` - Ratio of area at throat to area at throat. Always 1 (exists because f_ae also exists, f_ae being ratio of combustor area to throat area)
* `pip` - Ratio of pressure in chamber to pressure at exit
  * `t_ae` - Ratio of pressure in chamber to pressure at throat.


### Using ThermoInterface

Provided with this library is an interface to the thermo_spg.inp file provided with the library.

You can access materials as if `ThermoInterface` is a dictionary using materials' CEA Names. E.G. `ThermoInterface["CH4"]`. This object support checking for inclusion and iterations, such as `"CH4" in ThermoInterface` and `for name in ThermoInterface`. It supports dictionary methods such as `.keys()`, `.values()` and `.items()`.

The value returned by `ThermoInterface` accesses is a `ThermoMaterial object`, which is an `Output` object (dictionary that can also be accessed with .dot notation) with the following keys:

- `name` - Name of the material
- `reference` - Reference, if given. Otherwise ""
- `elements` - Dictionary of element: numerical value specified
- `condensed` - True if condensed phase, False otherwise
- `mw` - Molecular weight in g/mol
- `hf` - float heat of formation at 298.15 in kJ/mol (or assigned enthalpy if 0 temp range)
- `temp_ranges` - List of two-tuples of [range start, range end] (K)
- `reactant_only` - True if material only shows up in reactants, False otherwise

#### Available Methods

- `ThermoInterface.get_close_matches(name, [n])` - Gets close matches to a given material name. For example "Al(cr)" returns 'AL(cr)', 'ALN(cr)', 'Ag(cr)', 'W(cr)'. n is the maximum number of results returned (though fewer may be returned).
- `ThermoMaterial.defined_at(temp)` - Returns True if the material is specified at the given temperature, False otherwise. Materials specified at one temperature are actually allowed at that temperature +- 10K.

### Utilities

- ```open_thermo_lib()```: Opens the default thermo library input file using the user's default .inp file viewer (should prompt if none)

- ```open_pdfs()```: Opens the attached NASA pdfs using the user's default pdf viewer

- ```print_assets_directory()```: Prints to console the current location of the directory where CEA_Wrap assets are located. Also returns this value

- ```get_simple_thermo_lib_line(name, comment, atoms, is_condensed, mol_wt, hf, temperature=298.15)```: Returns and a string which can be inserted to represent a molecule in the Thermo Lib. Any entries which are longer than the allotted space results in an error
  
  - `name` (`str`): 24 chars max, Species name, such as CO or Fe2O3. These are assumed to be in a gas phase unless appended with (a) for agglomerate, (cr) for crystalline (solid), or (L) for liquid phase.
  - `comment` (`str`): 56 char max, Human-readable name and additional information
  - `atoms` (`dict[str, int]`): 5 atom max, A dictionary of "atomic symbol": number of atoms in molecule.

      Example: Water, "H2O", would be {"H": 2, "O": 1}

      Complex Example: Lead Acetate, "Pb(C2H3O2)2", would be {"Pb": 1, "C": 4, "H": 6, "O": 4}. Note that the C, H, and O are doubled because they are the *sum* of atoms in the molecule

      Note: "E" is a special value which represents an electron for ionic compounds. You can have negative amounts of "E" to indicate positively charged atoms

  - `is_condensed` (`bool`): True if the material is a solid/liquid, False if it is a gas.
    
  - `mol_wt` (`float`): Molecular weight of molecule in g/mol or kg/kmol
    
  - `hf` (`float`): Assigned Enthalpy of material at temperature specified. J/mol. If temperature=298.15K, this is the Heat of Formation.
    
  - `temperature` (`float`, *optional*, default=298.150): Temperature that the molecule is specified at, in Kelvin

- ```print_simple_thermo_lib_line(*args, **kwargs)```

  - A simple alias for `print(get_simple_thermo_lib_line(*args, **kwargs))`

- ```reload_thermo_lib()```: Moves the thermo and transport libs from your assets directory (if they have changed) and reloads ThermoInterface from the thermo_spg.inp in your assets directory
Use this if you have recompiled your thermo lib but do not want to restart python.

  Note: This is actually defined in CEA.py, not utils.py
  
#### DataCollector

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

### Environment Variables

CEA_Wrap offers a variety of environment variables to customize the locations and behavior when accessing necessary files. In the below options, a "Truthy" value is something like "true", "yes", or "1". A "Falsy" value is something like "false", "no", or "0". This determination only looks at the first character and is case insensitive.

Note: Most of these environment variables are only examined once on startup, so they must be set before the library is imported

- `CEA_USE_SITE_PACKAGES`: If set to some Truthy value, CEA_Wrap will not copy assets files to a local directory, instead using the site-packages where CEA_Wrap is installed. This may desired if using CEA_Wrap in a portable installation. This is not recommended if you make modifications to the thermo_spg.inp file, thermo.lib, or trans.lib, because these files will be overwritten if you update CEA_Wrap. DEFAULT: False
- `CEA_ASSETS_DIR`: If set, CEA_Wrap will use this directory as the "local" directory for finding "assets" such as the CEA executable, thermo.lib, and documentation PDFs. If this is set, the "appdirs" library dependency is not required. This is ignored if CEA_USE_SITE_PACKAGES is set. DEFAULT: None (uses appdirs.user_data_directory)
- `CEA_EXE_LOCATION`: If set, overrides the path to the CEA executable used for running CEA calculations
- `CEA_THERMO_LIB`: If set, overrides the path to the thermo.lib file used by the CEA executable
- `CEA_TRANS_LIB`: If set, overrides the path to the trans.lib file used by the CEA executable
- `CEA_THERMO_INP`: If set, overrides the path to the thermo_spg.inp file used by the ThermoInterface in CEA_Wrap
- `CEA_COPY_THERMO_TO_TEMP`: If set to some Truthy value, when a CEA object is instantiated (such as when you import CEA_Wrap), the thermo.lib and trans.lib from assets are copied to a temporary directory (by default). This enables multiprocessing because CEA holds a lock on the thermo.lib and trans.lib files while it is running, and will crash the program if two instances of CEA attempt to access the thermo.lib files at the same time. DEFAULT: True
- `CEA_USE_LEGACY`: Uses a modified CEA interface that uses legacy logic for interacting with the original FCEA2 executable

### Making A Portable Installation

To download the package to be used as a portable installation, simply clone/download this repository.

#### Installing/Using the portable installation

In whatever environment you are using, you can navigate to the folder containing `setup.py`, then call `pip install -e .`

This will install the current directory into your environment as an "editable" package. **Alternatively**, you can simply place your code in the same directory as the `setup.py` script. So your directory would look like:

```
My_Folder
│   my_code.py   <<<<
|   setup.py
|   setup.cfg
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

#### Setting the environment variable

Before you run your code, you must set the environment variable `CEA_USE_SITE_PACKAGES` to '1' or 'True' or similar. You must do this **before** you import CEA_Wrap. You can set the environment variable in a variety of ways. If your code is in "my_code.py", ...

##### Set environment variables in Windows:

```shell
C:\My_Folder> set CEA_USE_SITE_PACKAGES=1
C:\My_Folder> python my_code.py
```

##### Set environment variables in Linux:

```bash
~$ export CEA_USE_SITE_PACKAGES=1
~$ python3 my_code.py
```

##### Set environment variable from within code

`my_code.py`

```python
import os
os.environ["CEA_USE_SITE_PACKAGES"] = '1'
from CEA_Wrap import Oxidizer, Fuel, RocketProblem
...
```

### Adding materials/recompiling the thermo lib (advanced usage)
  
  You can add new molecules by modifying your thermo_spg.inp file and then recompiling the thermo lib.

  1. Find the location of your thermo_spg.inp file by using print_assets_directory()
  1. Open thermo_spg.inp in a text editor.
  1. Scroll to the bottom
  1. You can use print_simple_thermo_lib_line() to create a simple entry in the format of other reactants
  1. Add your entry below the other entries but above "END REACTANTS". Ensure your entry has similar spacing to existing entries
  1. To recompile your thermo lib, either use the "thermo_lib_recompile" batch file or run "FCEA2" and type "thermo_spg"
  1. If you have a python kernel open, you should now restart it or use the "reload_thermo_lib()" function to reload the thermo lib used.
