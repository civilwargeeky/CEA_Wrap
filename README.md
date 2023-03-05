# CEA_Wrap
A Python-Based wrapper for the NASA CEA Thermochemical Code

# Installation Instructions
We are now on [PyPi!](https://pypi.org/project/CEA-Wrap/)
In a command prompt type ```pip install --upgrade CEA_Wrap``` to upgrade/install CEA_Wrap

You can import it as any other python module with ```import CEA_Wrap```. Whenever you import the file, it will put the required thermo.lib and trans.lib files into your current directory.

## Installation on Mac
On mac, the installation script will attempt to compile the fortran executable on your system. Should it fail to do so, you will have to compile it manually.
You can see a discussion from a successful user here: [Issue #1](https://github.com/civilwargeeky/CEA_Wrap/issues/1#issuecomment-1033918162)

## Upgrading from Version <1.5.0
If you upgrade and your version is <1.5.0, any custom thermo lib you have made will be overwritten

Please make a copy of your thermo_spg.inp, place it in the new data directory, and recompile after upgrading

As of 1.5.0, thermo_spg.inp and all other assets are kept in a data directory, rather than the package assets directory. This means that your custom thermo lib will no longer be overwritten when you upgrade

# Examples
Examples on basic use can be found in the "examples" directory. Feel free to download them and try them out!

Or you can run a short demo by doing "python -m CEA_Wrap" on the command line!

# Documentation

## Materials
  In order to run problems, you must create materials. `Materials` must be either Fuel or Oxidizer (alias: F or O) objects
  
  ### Constructor Parameters - Either Fuel() or Oxidizer():
  ```Material(name, temp=298.15, wt_percent=None, mols=None, chemical_composition = None, hf = None):```
  
* `name`: required parameter, the CEA material name with correct spelling. E.G. aluminum is "AL(cr)" and methane is "CH4". 
  * If you specify chemical_composition, name can be whatever single word you want
  * If you try to specify a chemical which is not in the thermo_spg.inp file, a ValueError will be raised. To prevent this check, set `Material.check_against_thermo_inp = False`
* `temp`: default 298.15, specified reactant initial temperature, Kelvin
  * If you try to specify a temperature which is not supported for this material in the thermo_spg.inp file, a ValueError will be raised. To prevent this check, set `Material.check_against_thermo_inp = False`
* `wt_percent`: A weight-based percentage for the element. Weight percentages do not need to add up to 100, and are calculated on the ratio with other Fuels/Oxidizers
* `mols`: A mol-based percentage for the element. Can be used as in `Oxidizer("O2", mols=1)` and `Oxidizer("N2", mols=3.76)` for air (except CEA has "air" as a reactant...)
* **NOTE:** wt_percent and mols cannot be specified together, if neither is defined, the `Material` gets a wt_percent of 1
* `chemical_composition`: chemical composition such as "LI 1 B 1 H 4" for LiBH4. If defined, will not use CEA default values 
  * **NOTE: Not rigourously tested**
* `hf`:  Enthalpy of formation, kJ/mol, must be specified if chemical_composition is specified
  * **NOTE: Not rigourously tested**

### Available Members
* All parameters from the constructor are also members
* `.ref` - If Material.check_against_thermo_inp is True, this will be a ThermoMaterial representing the material

### Available Methods
* `.set_wt_percent(wt_percent)` - Sets the wt_percent for the `Material`. Sets mols to None
* `.set_mols(mols)` - Sets the mols for the `Material`. Sets wt_percent to None
* `.set_temp(temp)` - Sets the temp for the `Material`, in Kelvin. If you try to specify a temperature which is not support for this material, a ValueError will be raised
* `.is_mols()/.is_wt_percent()` - Checks if the `Material` is set in weight % or mol ratio

## Generic Problem Methods (apply to all problems such as DetonationProblem, RocketProblem, etc.):

### Constructor Parameters (making new Problem Objects):
  ```Problem(*, **kwargs)```:
* **NOTE:** all parameters must be specified by keyword, e.g. problem = RocketProblem(pressure=500, massf=True)
* `pressure`: default 1000, Initial problem pressure
* `materials`: default None, List of `Material` objects, order doesn't matter, of Oxidizer and Fuel objects e.g. materials=[material1, material2, ...]
  * materials can also be specified when you run a problem like `problem.run_cea(material1, material2, ...)`
  * materials MUST all have wt_percent specified or all have mols specified, can't have mixtures.
* `massf`: default False, CEA usually outputs product ratios in mole ratios, if massf is True, mass ratios are specified
* `filename`: default "my_output", the filename we save our .inp files to and CEA saves our .out and .plt to.
  * DO NOT INCLUDE ".inp" IN YOUR FILENAMES. 
* `pressure_units`: default "psi", the units that your input pressure is in. Possible values are "bar", "atm", "psi", or "mmh"
* `inserts`: default None, a list of CEA names for species which should be forced into the product considerations. Should be specified as either a space-separated string of names, a list of string names, or a list of `Material` objects.
  * Note: If you try to specify a chemical which is not in the thermo_spg.inp file, a ValueError will be raised. To prevent this check, set `Problem.check_against_thermo_inp =  False` 
  * Tip: If you are doing calculations with Aluminum, I recommend using inserts=["AL2O3(L)", "AL2O3(a)"]
* `omits`: default None, a list of CEA names for species which should be specifically ignored in the product considerations. Specified similar to inserts.
* **NOTE:** You must specify one of the following reactant ratio schemes before running a problem. Either during initialization with `x=1.0` or later with `problem.set_x(1.0)` where x is p_f, o_f, etc.
#### Specifying reactant ratios:
   Key  | CEA Key | Description
  ------|---------|---------------------------------------------------------------------------------------------------------
   p_f  |  %f     | Percent fuel by weight
   f_o  |  f/o    | Fuel-to-oxidant weight ratio
   o_f  |  o/f    | Oxidant-to-fuel weight ratio
   phi  |  phi    | Equivalence ratios in terms of fuel-to-oxidant weight ratios (eq. (9.19) in Gordon and McBride, 1994)
   r_eq |  r      | Chemical equivalence ratios in terms of valences (eq. (9.18) in Gordon and McBride, 1994)

### Available Methods
* `data = .run([*materials])`

Run the CEA problem, returning an "Output" object, which is similar to a dictionary (keys and values specified later in this documentation)
  * Inputs are optionally `Material`s to be used in this run. If materials are not specified as an initial parameter or with .set_materials, you can list them here.
  
* `.set_absolute_o_f()`

Have you ever specified all your components in absolute percentages, and then have to manually calculate the o_f ratio based on what is fuel and what is oxidizer?

Well no more! Just call this function and we calculate the correct o_f ratio for you so that your absolute percentages are correctly reflected in the problem.
Functions off of the existing material list, so call this after setting materials. This works by summing wt_percent for all oxidizers and then dividing by the same for fuels. Then it sets o/f to this value.

* `.set_pressure(pressure)` - sets pressure
* `.set_materials([material1, material2, ...])` - provided a list of materials, sets the materials list
* `.set_massf(massf)` - Sets massf to True or False
* `.set_inserts(inserts)` - Set inserts to a space-separated string or list of materials or `Material`s
* `.set_omits(omits)` - Sets omits, similar to inserts
* `.set_filename(filename)` - Sets problem filename
* `.set_pressure_units(units)` - Sets input pressure units
* `.set_p_f(pf)` - Sets % Fuel for problem
* `.set_f_o(f_o)` - Sets fuel/oxidizer ratio
* `.set_o_f(o_f)` - Sets oxidizer/fuel ratio
* `.set_phi(phi)` - Sets equivalence ratio
* `.set_r_eq(r_eq)` - Sets valence equivalence ratio

## Rocket Problem Constructor Additional Parameters:
For `RocketProblem(*, **kwargs)`
* `sup`: default 1, supersonic exit area/throat area ratio
  * sup can be specified later with .set_sup
* `sub`: default None, subsonic exit/throat area ratio
  * sub can be specified later with .set_sub
* **NOTE:** sup and sub cannot be specified at the same time
* `ae_at`: alias for `sup`
  * ae_at can be specified later with .set_ae_at
* `pip`: Pressure ratio of chamber pressure/exit pressure
  * pip can be specified later with .set_pip
* `fac_ma`: Finite Area Combustor, with mass flow (mdot) / combustor chamber area. Units of (kg/s)/m^2
  * fac_ma can be specified later with .set_fac_ma
  * Cannot be specified at the same time as fac_ac
* `fac_ac`: Finite Area Combustor, with ac/at: Ratio of combustor area to throat area
  * fac_ac can be specified later with .set_fac_ac
  * Cannot be specified at the same time as fac_ma
* `analysis_type`: default "equilibrium", whether to use equilibrium reactions or frozen. For using frozen specify "frozen" or "frozen nfz=1" for frozen at the chamber or "frozen nfz=2" for frozen at the throat
* `nfz`: default None, If `analysis_type` is "frozen", this will set the frozen location to the given point. 1 for chamber, 2 for throat
* `custom_nfz`: default None, If `analysis_type` is "frozen", this is the position within the nozzle that composition will be frozen at. Uses the same unit as ae/at or pip. Example: custom_nfz=2 to freeze composition at ae/at=2

### RocketProblem Methods
* `.set_sup(sup)` - Sets supersonic area ratio
* `.set_sub(sub)` - Sets subsonic area ratio
* `.set_ae_at(sup)` - Sets supersonic area ratio
* `.set_pip(pip)` - Sets pressure ratio
* `.set_analysis_type(analysis, nfz=None, custom_nfz=None)` - Sets analysis type, with optional frozen specifications as above
* `.set_fac_ma(fac)` - Sets finite area combustor, with mass flow/area ratio
* `.set_fac_ac(fac)` - Sets finite area combustor, with combustor/throat area ratio
* `.unset_fac()` - Unsets finite area combustor

## Available Output Dictionary Keys:
All Problem data objects are "Output" objects, which are similar to dictionaries, but can also be accessed with dot notation.

For example if you had "data = problem.run_cea()", and wanted pressure, you could do either data.p or data["p"]

In addition, all product dictionaries are also "Output" objects so to get H2O fraction, you could use data.prod_c.H2O or data["prod_c"]["H2O"] or data["prod_c"].H2O, etc.

### Detonation:
* `prod_c` - dictionary of chamber products, in mole or mass fractions (as specified in problem)
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
* `phi` - weight-based equivalence ratio of oxidizer/fuel
### HP (Specified Enthalpy and Pressure):
* `prod_c` - dictionary of chamber products, in mole or mass fractions (as specified in problem)
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
* `phi` - weight-based equivalence ratio of oxidizer/fuel
### Rocket:
* **NOTE : Properties are by default at exit. Chamber parameters are prefixed "c_" and throat properties "t_"**
* **NOTE : Properties not defined for frozen flow are marked with an asterisk (*)**
* **NOTE : All properties defined at the throat are also defined as "f_property" when Finite Area Combustor is enabled (defined fac_ac or fac_ma)**
* `prod_c` - dictionary of chamber products, in mole or mass fractions (as specified in problem)
* `*prod_t` - dictionary of throat products, in mole or mass fractions (as specified in problem)
* `*prod_e` - dictionary of exit products, in mole or mass fractions (as specified in problem)
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

## Using ThermoInterface
Provided with this library is an interface to the thermo_spg.inp file provided with the library. 

You can access materials as if `ThermoInterface` is a dictionary using materials' CEA Names. E.G. `ThermoInterface["CH4"]`. This object support checking for inclusion and iterations, such as `"CH4" in ThermoInterface` and `for name in ThermoInterface`. It supports dictionary methods such as `.keys()`, `.values()` and `.items()`.

The value returned by `ThermoInterface` accesses is a `ThermoMaterial object`, which is an `Output` object (dictionary that can be accessed with .dot notation) with the following keys:
* `name` - Name of the material
* `reference` - Reference, if given. Otherwise ""
* `elements` - Dictionary of element: numerical value specified
* `condensed` - True if condensed phase, False otherwise
* `mw` - Molecular weight in g/mol
* `hf` - float heat of formation at 298.15 in kJ/mol (or assigned enthalpy if 0 temp range)
* `temp_ranges` - List of two-tuples of [range start, range end] (K)
* `reactant_only` - True if material only shows up in reactants, False otherwise

### Available Methods
* `ThermoInterface.get_close_matches(name, [n])` - Gets close matches to a given material name. For example "Al(cr)" returns 'AL(cr)', 'ALN(cr)', 'Ag(cr)', 'W(cr)'. n influences the number of results returned, and is the maximum number of results returned.
* `ThermoMaterial.defined_at(temp)` - Returns True if the material is specified at the given temperature, False otherwise. Materials specified at one temperature are actually allowed at that temperature +- 10K.

## Utilities
  ```open_thermo_lib()```
  Opens the default thermo library input file using the user's default .inp file viewer (should prompt if none)
  
  ```open_pdfs()```
  Opens the attached NASA pdfs using the user's default pdf viewer
  
  ```print_assets_directory()```
  Prints to console the current location of the directory where CEA_Wrap assets are located. Also returns this value
  
  ```print_simple_thermo_lib_line(name, comment, atoms, isCond, molWt, hf)```
  Returns and prints a string which can be inserted to represent a molecule in the Thermo Lib. Any entries which are longer than the allotted space results in an error
  
  * `name`: 24 chars max, Species name, such as CO or Fe2O3. These are assumed to be in a gas phase unless appended with (a) for agglomerate, (cr) for crystalline (solid), or (L) for liquid phase.
  
  * `comment`: 56 char max, Human-readable name and additional information
  
  * `atoms`: 5 atom max, A dictionary of "atom symbol": number of atoms in molecule.
  
    Example: H2O would be {"H": 2, "O": 1}
    
    Special value is "E" which represents an electron for ionic compounds
    
  * `isCond`: True if condensed phase (non-gas), False otherwise
  
  * `molWt`: Molecular weight of molecule in g/mol or kg/kmol
  
  * `hf`: Heat of formation at 298.15 K, in J/mol
  
### DataCollector
  A DataCollector object will conveniently compile output from Problem outputs into list format. Example usage can be found in the "looping.py" example
  
  ```DataCollector(self, *args, keys=[], chamber_keys=[], throat_keys=[], exit_keys=[])```
  
* `keys`: Also accepts list of arguments, these are keys such as 'cond' or 't_cp' or 'c_p'
* `chamber_keys`: List of chamber species mol/mass fractions to be included in the output object. Ex: "H2O" or "CO2". The key in the ouptut will be the molecule name with "c_" prepended
* `throat_keys`: List of nozzle throat species mol/mass fractions to be included. The key in the output will be the molecule name with "t_" prepended
* `exit_keys`: List of nozzle exit species mol/mass fractions to be included. The key in the output will be the molecule name with nothing prepended.
  
  *Methods*
* `add_data(data)` - Data should be the output from Problem.run(). Appends the output of the keys specified in the initializer to the object
* `to_csv(filename, filename: str, keys:list=None, formatString="f"` - Writes the data to csv at **filename**, with the keys in **keys** (or all keys in this object if None)
     using **formatString** to format the csv entries
       
