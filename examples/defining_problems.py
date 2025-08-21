### This simple example file goes over the different ways to define problems, material ratios, and materials
from pprint import pprint

# Import only the things we'll need for this problem. You could also do from CEA_Wrap import * if you like polluting your global namespace
from CEA_Wrap import Fuel, Oxidizer, HPProblem, RocketProblem

#############################################################################
## First we'll show the different ways to define problems
#############################################################################

# For all these we will consider stoichiometric combustion of methane and oxygen
# Define our materials, weight percentages default to 100%
methane = Fuel("CH4")
o2 = Oxidizer("O2")

## NOTE: All these methods are equivalent and will produce identical output
## NOTE: The only things you need to specify for a problem are some form of o/f ratio and materials
##         all other parameters have defaults, such as a pressure of 1000 psi and a filename of "my_output"

## Problem Definition 1: All parameters in problem constructor
problem = HPProblem(phi=1, pressure=500, materials=[methane, o2], filename="HP_Problem")
data1 = problem.run()

## Problem Definition 2: Specify parameters with object
problem = HPProblem()
problem.set_phi(1)
problem.set_pressure(500)
problem.set_materials([methane, o2])
problem.set_filename("HP_Problem")
data2 = problem.run()

## Problem Definition 3: Mix and match. Any parameter can be set in constructor or from object
##   This is particularly useful when you are modifying a problem within a loop, such as looping through pressures
problem = HPProblem(materials=[methane, o2], filename="HP_Problem")
problem.set_phi(1)
problem.set_pressure(500)
data3 = problem.run()

## Problem Definition 4: As the others, but you can also specify materials when you call the problem
problem = HPProblem(phi=1, pressure=500, filename="HP_Problem")
data4 = problem.run(methane, o2)

# Let's double check they are all the same!
assert data1 == data2 == data3 == data4

# Now we'll print out the data so you can see what it produces when you run it
print("Displaying results for stoichiometric methane/oxygen combustion")
pprint(data1)

#############################################################################
### Now we'll show different ways to define material ratios and temperatures
#############################################################################

## NOTE: percentages/mol ratios are intra-fuels and intra-oxidizers. Ratios between fuels and oxidizers are covered in the next section

## Material Example 1: Defining by weight percentage
methane = Fuel("CH4", wt=20)
octane = Fuel("C8H18,isooctane", wt=80)
oxygen = Oxidizer("O2") # defaults to wt=100

data1 = HPProblem(phi=1).run(methane, octane, oxygen) # Stoichiometric between fuels and oxidizers, defaults 1000 psi

## Material Example 2: Defining by mol ratio
## NOTE: In CEA you can also specify "Air" as an oxidizer, so this is a poor way to specify this problem for example's sake
oxygen = Oxidizer("O2", mols=1)
nitrogen = Oxidizer("N2", mols=3.76)
methane = Fuel("CH4", mols=1) # Note: all materials must use mols or wt, no mixing, so we must specify mols=1
data2 = HPProblem(phi=1).run(oxygen, nitrogen, methane) # Stoichiometric methane/

## Material Example 3: Changing ratios after definition. This is especially useful in loops for materials
methane = Fuel("CH4", wt=20)
octane = Fuel("C8H18,isooctane", wt=80)
oxygen = Oxidizer("O2") # defaults to wt=100

# Define our problem
problem = HPProblem(phi=1, materials=[methane, octane, oxygen])
data3_1 = problem.run() # Run with original ratio

# Change to a 50/50 ratio of the two
## NOTE: percentages do not need to add up to 100%! CEA calculates from sum of fuel/oxidizer percentages
methane.wt = 1
octane.wt = 1
data3_2 = problem.run() # Run with new updated ratio

## Material Example 4: Specifying with temperatures
methane = Fuel("CH4")
o2 = Oxidizer("O2", temp=200) # Cool our oxygen a bit
problem = HPProblem(phi=1, materials=[methane, o2])
data4_1 = problem.run() # Run with this temp

o2.set_temp(5000) # then really hot oxygen
data4_2 = problem.run() # run with this temp

#############################################################################
## Lastly we'll show the different ways to specify to oxidizer-fuel ratios
#############################################################################

## Ratios Example 1: Equivalence ratios
air = Oxidizer("Air")
methane = Fuel("CH4")

# Define the problem with fuel-lean combustion
problem = HPProblem(materials=[methane, air], phi=0.7)
data1_1 = problem.run()
# Then we'll set fuel-rich
problem.set_phi(1.5)
data1_2 = problem.run()

## Ratios Example 2: O/F ratio
# Define the problem with slightly fuel-rich combustion (using materials from previous)
problem = HPProblem(materials=[methane, air], o_f=20)
data2_1 = problem.run()
# Then we'll set fuel-lean
problem.set_o_f(10)
data2_2 = problem.run()

## Ratios Example 3: Helper function: absolute o/f ratios
# So doing rocket problems with absolute percentages specified (not as percentage of fuel/oxidizer and o/f ratio specified)
#   I got annoyed, so I wrote a helper function that calculates the o/f ratio from wt of all constituents
# So we'll be doing a rocket problem with AP/Aluminum/HTPB
# Define materials
aluminum = Fuel("AL(cr)", wt=12) # (cr) for "crystalline" or condensed phase
htpb = Fuel("HTPB", wt=14) # This was added at Purdue so doesn't include (cr) in the name
ap = Oxidizer("NH4CLO4(I)", wt=74) # ammonium perchlorate (form I, specified at room temperature)

# Specify a RocketProblem with supersonic expansion ratio of 15 at chamber pressure of 1500 psi
# We don't specify our phi or o_f as we want it to calculate that for us
problem = RocketProblem(pressure=1500, sup=15, materials=[aluminum, htpb, ap])
# sets the o_f ratio by summing wt for each oxidizer and  fuel
#  and dividing by the sum for fuels
# If you change material wts, remember to call this again!
problem.set_absolute_o_f() # sets o/f to (74)/(12+14)
data4 = problem.run()

## Ratios Example 4: miscellaneous other methods
# Specify problem by percent fuel (by weight)
problem = HPProblem(materials=[methane, air], p_f=20)
problem.set_p_f(30)
# Specify problem by f/o ratio (because sure, why not)
problem = HPProblem(materials=[methane, air], f_o=1/20) # fuel-rich
problem.set_f_o(1/10) # fuel-lean
# Specify problem by "Chemical equivalence ratios in terms of valences" (I don't actually know what this is)
problem = HPProblem(materials=[methane, air], r_eq=1) # stoichiometric by valences, I guess
problem.set_r_eq(1.5) # Set it to something else
