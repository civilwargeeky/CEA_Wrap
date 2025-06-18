# Example for how to loop through materials and collate the properties you are concerned about
# Note: To run this example you will need numpy and matplotlib

try:
  import numpy as np
  import matplotlib.pyplot as plt
except ImportError:
  print("You will need to install matplotlib and numpy to run this example")
  import sys
  sys.exit(1)
from time import time

# Only import things we'll need for this problem
from CEA_Wrap import Fuel, Oxidizer, RocketProblem, DataCollector

# For all the examples in this problem we'll use Aluminized AP/HTPB as the reactants
aluminum = Fuel("AL(cr)", wt=12) # (cr) for "crystalline" or condensed phase
htpb = Fuel("HTPB", wt=14) # This was added at Purdue so doesn't include (cr) in the name
ap = Oxidizer("NH4CLO4(I)", wt=74) # ammonium perchlorate (form I, specified at room temperature)
m_list = [aluminum, htpb, ap] # for convenience so I can pass it into all problems

# We'll actually use the same problem for the first 3 examples
problem = RocketProblem(materials=m_list)
problem.set_absolute_o_f() # have it calculate o_f for us from material percentage

## Example 1: Using DataCollector objects (Recommended)
print("Running Data Collector Case")
pressures = np.linspace(1000, 2000, 25) # pressures from 1000 to 2000
collector = DataCollector("c_t") # want to collect chamber temperature
for pressure in pressures:
  problem.set_pressure(pressure)
  collector.add_data(problem.run()) # adds all desired variables to collector
  
plt.plot(pressures, collector.c_t) # could also do collector["c_t"] for the same effect
plt.title("Collector Pressure vs Temperature")
plt.xlabel("Pressure (psi)")
plt.ylabel("Chamber Temperature (K)")
plt.show()


## Example 2: Using standard python arrays
print("Running Python Array Case")
pressures = [1000 + i*40 for i in range(25)] # pressures from 1000 to 2000
temperatures = []
for pressure in pressures:
  problem.set_pressure(pressure)
  data = problem.run()
  temperatures.append(data.c_t) # append chamber temperature

plt.plot(pressures, temperatures)
plt.title("Python Arrays Pressure vs Temperature")
plt.xlabel("Pressure (psi)")
plt.ylabel("Chamber Temperature (K)")
plt.show()


## Example 3: Using numpy arrays
print("Running Numpy Array Case")
arr_size = 25 # 25 runs of CEA
pressures = np.linspace(1000, 2000, arr_size) # pressures from 1000 to 2000
temperatures = np.zeros(arr_size) # initialize with 0s
for i, pressure in enumerate(pressures):
  problem.set_pressure(pressure)
  data = problem.run()
  temperatures[i] = data.c_t # set current chamber temperature

plt.plot(pressures, temperatures)
plt.title("Numpy Pressure vs Temperature")
plt.xlabel("Pressure (psi)")
plt.ylabel("Chamber Temperature (K)")
plt.show()


## Example 4: Terrible no good using numpy arrays to get all your properties
print("Running terrible Numpy Array Case")
arr_size = 30 # runs of CEA
pressures = np.linspace(1000, 2000, arr_size) # pressures from 1000 to 2000
# Define our empty arrays for every property
# It would be the same for python arrays
t = np.zeros(arr_size)
isp = np.zeros(arr_size)
cf = np.zeros(arr_size)
gamma = np.zeros(arr_size)
cstar = np.zeros(arr_size)
mach = np.zeros(arr_size)
exit_cp = np.zeros(arr_size)
chamber_cp = np.zeros(arr_size)
percent_hcl = np.zeros(arr_size)
percent_alumina = np.zeros(arr_size)
# Run the problem at different pressures
for i, pressure in enumerate(pressures):
  problem.set_pressure(pressure)
  data = problem.run()
  t[i] = data.c_t # chamber temperature
  isp[i] = data.ivac # vacuum isp
  cf[i] = data.cf # exit thrust coefficient
  gamma[i] = data.gamma # exit real ratio of specific heats
  cstar[i] = data.cstar # characteristic velocity
  mach[i] = data.mach # exit mach number
  exit_cp[i] = data.cp # exit specific heat
  chamber_cp[i] = data.c_cp # chamber specific heat
  percent_hcl[i] = data.prod_e.HCL
  try: # Fun fact: there may not be any liquid alumina in the exhaust products, so you have to do this
       #    For every single product you want to look at.
       #    blech
    percent_alumina[i] = data.prod_e["AL2O3(L)"]
  except KeyError:
    percent_alumina[i] = 0
# I'm not going to plot these, but you get the idea.


## Example 5: Using the awesome DataCollector objects to make your life easy
print("Running the awesome super cool Data Collector Case")
pressures = np.linspace(1000, 2000, 30) # pressures from 1000 to 2000
# Define a collector to collect all these properties. Also get composition at the exit for these product species
#   If you wanted to define keys for the chamber, use "chamber_keys" instead
#   NOTE: You cannot have chamber and exit keys in the same collector. If you want both, use multiple collectors
collector = DataCollector("c_t", "ivac", "cf", "gamma", "cstar", "mach", "cp", "c_cp", exit_keys=["HCL", "AL2O3(L)"])
for pressure in pressures:
  problem.set_pressure(pressure)
  collector.add_data(problem.run()) # adds all desired variables to collector
# And that's it! Now we can access our array of chamber temperatures with
collector.c_t
# And our mol_ratios of HCL in the exhaust with
collector.HCL
# The really cool thing is that we don't need to worry about species not existing, because they will automatically be
#   set to 0 in the correct position of the array
# Fun thing! You can also write your collectors to csv
collector.to_csv("looper.csv")


## Example 6: Varying ratios of materials
print("Running various aluminum contents (may take a bit)")
percent_aluminum = np.linspace(0, 25, 250) # run percentage of aluminum from 0% to 15%
# reminder: 12% HTPB, and we'll reduce the amount of AP to compensate
collector = DataCollector("c_t", "ivac") # Show chamber temperature and isp dependence on % aluminum

start_time = time()
for al_percent in percent_aluminum:
  aluminum.wt = al_percent
  # NOTE: CEA doesn't like it when a material has 0%, so at 0% the material is not entered into the .inp file
  ap.wt = (100-12-al_percent)
  problem.set_absolute_o_f() # change our o/f to reflect situation
  collector.add_data(problem.run())
end_time = time()
print("Completed in {}s! {} CEA calls per second".format(end_time-start_time, 250/(end_time-start_time)))

# Adapted from https://matplotlib.org/stable/gallery/subplots_axes_and_figures/two_scales.html
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Percentage Aluminum (%)')
ax1.set_ylabel('Temperature', color=color)
ax1.plot(percent_aluminum, collector.c_t, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('Vacuum Isp (s)', color=color)
ax2.plot(percent_aluminum, collector.ivac, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# Example 7: Trying the same thing threaded
print("Running super-fast multi-threaded case")
from concurrent.futures import ThreadPoolExecutor

collector = DataCollector("c_t", "ivac") # Show chamber temperature and isp dependence on % aluminum

def inner_function(al_percent): # This function is just the inside of the previous loop
  # It is wise to re-define our materials and function for each call to this function, so that different runs do not over-write each other's material amounts
  aluminum = Fuel("AL(cr)", wt=12) # (cr) for "crystalline" or condensed phase
  htpb = Fuel("HTPB", wt=14) # This was added at Purdue so doesn't include (cr) in the name
  ap = Oxidizer("NH4CLO4(I)", wt=74) # ammonium perchlorate (form I, specified at room temperature)
  aluminum.wt = al_percent
  # NOTE: CEA doesn't like it when a material has 0%, so at 0% the material is not entered into the .inp file
  ap.wt = (100-12-al_percent)
  # We also need to re-define our problem in each function so that we do not accidentally over-write problem components in other threads
  problem = RocketProblem(materials=[aluminum, htpb, ap])
  problem.set_absolute_o_f() # change our o/f to reflect situation
  collector.add_data(problem.run(), al_percent) # Because the threaded cases may run out of order, we add our independent variable to our data collector, so we can sort later

start_time = time()
with ThreadPoolExecutor() as thread_pool:
  thread_pool.map(inner_function, percent_aluminum)
collector.sort() # Sort our data collector by our independent variable
end_time = time()
print("Completed in {}s! {} CEA calls per second".format(end_time-start_time, 250/(end_time-start_time)))


# # Adapted from https://matplotlib.org/stable/gallery/subplots_axes_and_figures/two_scales.html
# fig, ax1 = plt.subplots()
# color = 'tab:red'
# ax1.set_xlabel('Percentage Aluminum (%)')
# ax1.set_ylabel('Temperature', color=color)
# ax1.plot(percent_aluminum, collector.c_t, color=color)
# ax1.tick_params(axis='y', labelcolor=color)
# ax2 = ax1.twinx()
# color = 'tab:blue'
# ax2.set_ylabel('Vacuum Isp (s)', color=color)
# ax2.plot(percent_aluminum, collector.ivac, color=color)
# ax2.tick_params(axis='y', labelcolor=color)
# fig.tight_layout()  # otherwise the right y-label is slightly clipped
# plt.show()