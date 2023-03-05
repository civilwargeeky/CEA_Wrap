# Import the things we need to do our task. "F" and "O" are aliases for fuel and oxidizer
from CEA_Wrap import Fuel, Oxidizer, RocketProblem

#############################################################################
## Example 1: Simple 1 fuel 1 oxidizer liquid rocket
#############################################################################

# Cyrogenic fuels
h2 = Fuel("H2(L)", temp=20) # Liquid H2, Temp in Kelvin
lox = Oxidizer("O2(L)", temp=90)
# Rocket at 2000psi and supersonic area ratio of 5. Frozen in chamber
problem1 = RocketProblem(pressure=2000, materials=[h2, lox], phi=1, sup=5, analysis_type="frozen nfz=1")
results = problem1.run()

# For a full listing of available members, see the documentation at https://github.com/civilwargeeky/CEA_Wrap
print("Stoichiometric, cryogenic rocket. Frozen at chamber")
print("Temperature at exit (K):", results.t)
print("Ratio of Specific Heats:", results.gamma)
print("Vacuum Isp:             ", results.ivac)

# Now let's freeze it at the throat instead
problem1.set_analysis_type("frozen", nfz=2)
results = problem1.run()

print("Stoichiometric, cryogenic rocket. Frozen at throat")
print("Temperature at exit (K):", results.t)
print("Ratio of Specific Heats:", results.gamma)
print("Vacuum Isp:             ", results.ivac)


# Now let's freeze it at a custom point within the nozzle
# This uses the same units as our specified exit condition, so could be using sup, sub, or pip
problem1.set_analysis_type("frozen", custom_nfz=1.5) # Freeze it at ae/at == 1.5
results = problem1.run()

print("Stoichiometric, cryogenic rocket. Frozen at ae/at = 1.5")
print("Temperature at exit (K):", results.t)
print("Ratio of Specific Heats:", results.gamma)
print("Vacuum Isp:             ", results.ivac)

# Comparison to Equilibrium
problem1.set_analysis_type("equilibrium")
results = problem1.run()

print("Stoichiometric, cryogenic rocket. Equilibrium composition")
print("Temperature at exit (K):", results.t)
print("Ratio of Specific Heats:", results.gamma)
print("Vacuum Isp:             ", results.ivac)