from .CEA import Fuel, Oxidizer, HPProblem, RocketProblem
import pprint

code1 = """
# Define our materials, weight percentages default to 100%
methane = Fuel("CH4")
o2 = Oxidizer("O2")
# Define our problem with stoichiometric combustion
# Default pressure unit is psi, but can be specified with pressure_unit=
problem = HPProblem(phi=1, pressure=500, materials=[methane, o2], filename="HP_Problem")
data = problem.run() # Actually run CEA
"""

print("We will first run a simple HP Problem (specified enthalpies and pressure)")
print("The code used is:\n")
print("-"*79)
print(code1)
print("-"*79,"\n")
print("Press Enter to Continue...")
input()

# Define our materials, weight percentages default to 100%
methane = Fuel("CH4")
o2 = Oxidizer("O2")
# Define our problem with stoichiometric combustion
# Default pressure unit is psi, but can be specified with pressure_unit=
problem = HPProblem(phi=1, pressure=500, materials=[methane, o2], filename="HP_Problem")
data = problem.run() # Actually run CEA

print("Which returns:")
pprint.pprint(data)

print("Press Enter to Continue...")
input()
print("-"*79,"\n")
print("-"*79,"\n")

code2 = """
# Define materials
# at 400K, 12% of total composition
aluminum = Fuel("AL(cr)", wt_percent=12, temp=400)
# htpb is only specified in our thermo lib at 298K +- 15K
htpb = Fuel("HTPB", wt_percent=14)
# ammonium perchlorate
ap = Oxidizer("NH4CLO4(I)", wt_percent=74, temp=400)

# Specify a RocketProblem with supersonic expansion ratio of 15 
#   at chamber pressure of 1500 psi
# We don't specify our phi or o_f as we want it to calculate that for us
problem = RocketProblem(pressure=1500, sup=15, materials=[aluminum, htpb, ap])
# sets the o_f ratio by summing wt_percent for each oxidizer and  fuel
#  and dividing by the sum for fuels
# If you change material wt_percents, remember to call this again!
problem.set_absolute_o_f()

data = problem.run()
"""

print("Next we will run a rocket problem with Aluminized AP/HTPB\n  at an elevated temperature")
print("The code used is:\n")
print("-"*79)
print(code2)
print("-"*79)
print("Press Enter to Continue...")
input()

# at 400K, 12% of total composition
aluminum = Fuel("AL(cr)", wt_percent=12, temp=400)
# htpb is only specified in our thermo lib at 298K +- 15K
htpb = Fuel("HTPB", wt_percent=14)
# ammonium perchlorate
ap = Oxidizer("NH4CLO4(I)", wt_percent=74, temp=400)
# Specify a RocketProblem with supersonic expansion ratio of 15 
#   at chamber pressure of 1500 psi
# We don't specify our phi or o_f as we want it to calculate that for us
problem = RocketProblem(pressure=1500, sup=15, materials=[aluminum, htpb, ap])
# sets the o_f ratio by summing wt_percent for each oxidizer and  fuel
#  and dividing by the sum for fuels
# If you change material wt_percents, remember to call this again!
problem.set_absolute_o_f() 
data = problem.run()

print("Which gives us:")
pprint.pprint(data)