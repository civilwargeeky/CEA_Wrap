### This example file goes over using the supplied thermo interface functionality

from CEA_Wrap import Fuel, Oxidizer, HPProblem, ThermoInterface

## Example 1: ThermoInterface saving you from making input errors
try:
  ap = Oxidizer("NH4CL04") # made two mistakes, used "0" instead of "O" and didn't specify form (I) or (II)
except ValueError as e:
  print("Error!", e)
  
# You can also turn this behavior off for whatever reason
Oxidizer.check_against_thermo_inp = False
Oxidizer("NH4CL04") # Now it allows you to do this 

Oxidizer.check_against_thermo_inp = True
# NOTE: This also works for inserts and omits on Problems, and can be turned off in the same way

## Example 2: Accessing properties of output species
methane = Fuel("CH4")
air = Oxidizer("Air")
data = HPProblem(pressure=14.7, materials=[methane, air], phi=1).run() # Create and run problem
print(f'{"Material":^10}|{"Mol Frac":^10}|{"MW (g/mol)":^10}|{"Hf (kJ/mol)":^12}|{"Condensed?":^10}') # Headings spaced out
for element in sorted(data.prod_c): # Goes through the keys of output species
  dat = ThermoInterface[element]
  print("{:>10}|{:^10.4f}|{:10.3f}|{:12.3f}|{}".format(element, data.prod_c[element], dat.mw, dat.hf, dat.condensed))