# Import the things we need to do our task. "F" and "O" are aliases for fuel and oxidizer
from CEA_Wrap import Fuel, Oxidizer, RocketProblem

#############################################################################
## Example 1: Simple 1 fuel 1 oxidizer liquid rocket
#############################################################################

# Cyrogenic fuels
h2 = Fuel("H2(L)", temp=20) # Liquid H2, Temp in Kelvin
lox = Oxidizer("O2(L)", temp=90)
# Rocket at 2000psi and supersonic area ratio of 5
problem1 = RocketProblem(pressure=2000, materials=[h2, lox], phi=1, sup=5)
results = problem1.run()

# For a full listing of available members, see the documentation at https://github.com/civilwargeeky/CEA_Wrap
print("Stoichiometric, cryogenic rocket")
exit_pressure    = results.p
chamber_pressure = results.c_p
exit_cp          = results.cp
chamber_cp       = results.c_cp
exit_isp         = results.isp
throat_isp       = results.t_isp
print("Pressures (bar):", exit_pressure, chamber_pressure)
print("Cp (kJ/(kg*K):", exit_cp, chamber_cp)
print("Isp (s):", exit_isp, throat_isp)

# We can also access the exhaust products (in default mol fraction)
percent_water = results.prod_e.H2O # Can also do results.prod_e["H2O"]
print("Percent Water in exhaust:", percent_water)

# Now, let's increase our propellant temperature and see how that impacts isp
# Because we are going above the material boiling points, we will have to switch to gaseous hydrogen and oxygen
h2 = Fuel("H2") # Default temperature is 297
o2 = Oxidizer("O2")
problem1.set_materials([h2, o2])
results = problem1.run()
print("Hotter Isp:", results.isp, results.t_isp)

# Now let's increase our propellant temperature even more!
# We can just modify the materials in-place since they are already gas phase
h2.set_temp(800) # 500C
o2.set_temp(800)
results = problem1.run()
print("Even Hotter Isp:", results.isp, results.t_isp)
print("Wow, that didn't change very much...")

#############################################################################
## Example 2: Multiple fuel/oxidizer solid rocket
#############################################################################

# Let's make a solid rocket with multiple fuels and oxidizers
# It's composition will be 10% Aluminum, 50% AP, 20% AN, and 20% HTPB
print("\n\nDoing Solid Rocket")

### METHOD 1: Calculate the o/f yourself
AL = Fuel("AL(cr)", wt_percent=1/3) # AL is 1/3 of the fuels
HTPB = Fuel("HTPB", wt_percent=2/3) # HTPB is 2/3 of the fuels
AP = Oxidizer("NH4CLO4(I)", wt_percent=50/(50+20)) # The I and IV are different entries of the same molecule specified at different temperatures
AN = Oxidizer("NH4NO3(IV)", wt_percent=20/(50+20))

problem2 = RocketProblem(pressure=1000, materials=[AL, HTPB, AP, AN], o_f = (50+20)/(10+20), sup=5)
results = problem2.run()
print("Temperature (K), Isp (s) =", results.t, results.isp)

### METHOD 2: Automatically Calculate o_f
AL = Fuel("AL(cr)", wt_percent=10) # 10% of total
HTPB = Fuel("HTPB", wt_percent=20) 
AP = Oxidizer("NH4CLO4(I)", wt_percent=50) # The I and IV are different entries of the same molecule specified at different temperatures
AN = Oxidizer("NH4NO3(IV)", wt_percent=20)

problem2 = RocketProblem(pressure=1000, materials=[AL, HTPB, AP, AN], sup=5)
problem2.set_absolute_o_f() # Calculates o_f assuming that material wt_percents are relative to full mixture
results = problem2.run()
print("Temperature (K), Isp (s) =", results.t, results.isp)

# We can also examine the exhaust products (sorted by value instead of key)
product_tuples = results.prod_e.items()
for key, value in sorted(product_tuples, key=lambda value: value[1], reverse=True):
  print("{:20} : {:8.4%}".format(key, value))