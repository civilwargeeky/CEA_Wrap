# Example for how to loop through materials and collate the properties you are concerned about
# Note: To run this example you will need numpy and matplotlib

try:
  import numpy as np
  from scipy.optimize import minimize, LinearConstraint, approx_fprime
except ImportError:
  print("You will need to install scipy and numpy to run this example")
  import sys
  sys.exit(1)

# Only import things we'll need for this problem
from CEA_Wrap import Fuel, Oxidizer, RocketProblem

# We will set a few constants up here so you can easily play with different numbers
CONST_WT_HTPB = 14 # Constant for our problem-1 propellant and initial guess in problem 2
EXIT_RATIO = 2 # Ratio of exhaust area to throat area
PRESSURE = 1000 # psi

# Optimizer-specific options
INITIAL_AL = 1 # % Aluminum to guess
EPSILON = 0.2 # In weight %: Step size used in the numerical approximation of the jacobian

# For all the examples in this script we'll use Aluminized AP/HTPB as the reactants
aluminum = Fuel("AL(cr)") # (cr) for "crystalline" or condensed phase
ap = Oxidizer("NH4CLO4(I)") # ammonium perchlorate (form I, specified at room temperature)
htpb = Fuel("HTPB", wt=CONST_WT_HTPB) # This was added at Purdue so doesn't include (cr) in the name

problem = RocketProblem(pressure=PRESSURE, materials=[aluminum, htpb, ap], sup=EXIT_RATIO,
                        phi=1) # We set an arbitrary fuel-oxidizer ratio, as we will need to change it later

### Problem Statement: We want to find a propellant that will *maximize* Isp in an Aluminum + AP + HTPB propellant.
### Solution 1: An example of a simple "optimization" that may get you through many problems. 

aluminum_amounts = np.arange(0, 30, step=0.5) # We suspect that the optimum should be between 0 and 30% aluminum
isp_outs = np.zeros_like(aluminum_amounts) # Set up our array of isp values

for i, wt_aluminum in enumerate(aluminum_amounts):
  aluminum.set_wt_percent(wt_aluminum) # Update our material
  ap.set_wt_percent(100 - CONST_WT_HTPB - wt_aluminum) # Update our AP to make our propellant sum to 100% propellant
  problem.set_absolute_o_f() # Set the correct o_f ratio
  results = problem.run() # Re-run the problem

  isp_outs[i] = results.isp

ideal_index = np.argmax(isp_outs) # index of maximum Isp
ideal_isp = isp_outs[ideal_index]
ideal_al = aluminum_amounts[ideal_index]
print("-- Case: Simple Enumeration of Options -- ")
print(f"Optimal Amount of Aluminum: {ideal_al}")
print(f"Ideal Isp at that weight:   {ideal_isp}")


### Solution 2: We will modify the amount of aluminum, keep the HTPB constant, and change the AP to match 100
initial_values = [INITIAL_AL] # Starting wt percent guess of aluminum
bounds = [(0, 100-CONST_WT_HTPB)] # Set bounds for wt percent of aluminum. More than (100-74) would make negative AP. Note that it is a `list` of `tuples`

# Define our optimization function. This will get a `tuple` of aluminum weight percent set by the optimizer. The tuple is of size 1.
def to_optimize(x_tuple):
  wt_aluminum, = x_tuple # Putting the comma after `wt_aluminum` means "unpack the first element of the tuple into the variable wt_aluminum. Discard other values"
  # print("Examining", wt_aluminum) # Uncomment this to see what weights are explored
  aluminum.set_wt_percent(wt_aluminum) # Update our material
  ap.set_wt_percent(100 - CONST_WT_HTPB - wt_aluminum) # Update our AP to make our propellant sum to 100% propellant
  problem.set_absolute_o_f() # Set the correct o_f ratio
  results = problem.run() # Re-run the problem
  
  # Now, we set a "score" criteria, which the algorithm will try to *minimize*
  # As we want to *maximize* Isp, we will take the inverse
  return -results.isp # Note: We return *1* floating point number regardless of how many variables we have

optimization_result = minimize( # Documentation: https://docs.scipy.org/doc/scipy/reference/optimize.html
  to_optimize, # Function to minimize
  initial_values, # Initial values for each variable
  bounds=bounds, # Bounds for each variable
  options=dict( # Options provided to the solver. Our solver is "L-BFGS-B" because we have bounds but not constraints.
    # Options are here: https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html
    # Because the numbers from CEA are insensitive to tiny changes in weight %, we can't use the default forward-difference step size of 1e-8.
    # We set this to a larger number so it uses larger changes to examine the derivative in our scoring function
    eps=EPSILON
  )
)

ideal_al = optimization_result.x[0] # The `.x` is a list of variables in the optimal solution. We only have [0]: Wt Aluminum
ideal_isp = -optimization_result.fun # Value of the optimization function at the optimal solution. Invert it to get Isp instead of score
print("-- Case: 1D Optimization with Bounds -- ")
print(f"Optimal Amount of Aluminum: {ideal_al}")
print(f"Ideal Isp at that weight:   {ideal_isp}")


### Solution 3: Multi-variate optimization! Let's optimize every part of our propellant to see what works the best
# We will optimize both the Aluminum and HTPB, and calculate the AP from that

initial_values = [INITIAL_AL, CONST_WT_HTPB] # Starting wt percent guess of aluminum and HTPB
bounds = [ # Set bounds on reasonable amounts of all our constituents
  (0, 70), # For weight % of aluminum. We will just arbitrarily say 70% is the max
  (0, 70), # For weight % of HTPB. Again, arbitrary
]

# For this multi-variate problem, we *may* want to set up a "constraint" that our propellant should not have more than 100% propellant
# There are several ways to do this, I will use a LinearConstraint object. Documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.LinearConstraint.html
# We have 1 constraint that "Wt Aluminum" + "Wt HTPB" <= 100. We don't set a lower bound on the constraint as we can't have less than 0% of either
# So our constraint says "1 * Wt Aluminum" + "1 * Wt HTPB" <= 100
constraint = LinearConstraint(A=[[1,1]], ub=100)

# Define our optimization function. This will get a `tuple` of aluminum weight percent and htpb weight percent set by the optimizer. The tuple is of size 2.
def to_optimize(x_tuple):
  wt_aluminum, wt_htpb = x_tuple # Because x_tuple now has two elements, we unpack it into two variables
  # print(f"Examining {100-wt_aluminum-wt_htpb:.4f}% AP, {wt_aluminum:.4f}% Aluminum, and {wt_htpb:.4f}% HTPB") # Uncomment this to see what weights are explored
  aluminum.set_wt_percent(wt_aluminum) # Update our material
  htpb.set_wt_percent(wt_htpb)
  wt_ap = max(0, 100 - wt_htpb - wt_aluminum) # It is possible that AL and HTPB sum to >100 before our constraint gets applied. So we set AP to 0 in those cases
  ap.set_wt_percent(wt_ap) # Update our AP to make our propellant sum to 100% propellant
  problem.set_absolute_o_f() # Set the correct o_f ratio
  try:
    results = problem.run() # Re-run the problem
  except RuntimeError:
    # Sometimes, if you have too much aluminum, the nozzle will "freeze" and CEA fails to converge.
    # In this case, we still have to return something, so we will do a trick: We just return the positive amount of aluminum, so that the
    #   gradient of error-producing runs will prioritize less aluminum
    return wt_aluminum
  
  # Now, we set a "score" criteria, which the algorithm will try to *minimize*
  # As we want to *maximize* Isp, we will take the inverse
  return -results.isp # Note: We return *1* floating point number regardless of how many variables we have

# Our solver is "trust-constr" because we have bounds and `LinearConstraint`-based constraints not constraints.
# Very annoyingly, the trust-constr solver we does not implement "eps" to set the % change of components to calculate the derivative
#   Therefore, we must provide our own Jacobian with a specified difference
def approx_jacobian(x_tuple):
  return approx_fprime(x_tuple, to_optimize, epsilon=EPSILON)

optimization_result = minimize( # Documentation: https://docs.scipy.org/doc/scipy/reference/optimize.html
  to_optimize, # Function to minimize
  initial_values, # Initial values for each variable
  bounds=bounds, # Bounds for each variable
  constraints=[constraint], # List of constraints (we only have 1)
  jac=approx_jacobian # Our jacobian function, defined earlier
)

ideal_al = optimization_result.x[0] # The `.x` is a list of variables in the optimal solution. [0]: Wt Aluminum
ideal_htpb = optimization_result.x[1] # [1] is Wt HTPB
ideal_isp = -optimization_result.fun # Value of the optimization function at the optimal solution. Invert it to get Isp instead of score
print("-- Case: 2D Optimization with Bounds and Constraints -- ")
print(f"Optimal Amount of Aluminum: {ideal_al}")
print(f"Optimal Amount of HTPB:     {ideal_htpb}")
print(f"Optimal Amount of AP:       {100-ideal_al-ideal_htpb}")
print(f"Ideal Isp at that weight:   {ideal_isp}")