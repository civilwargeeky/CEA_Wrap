from CEA_Wrap import Rocket_Problem, Data_Collector, Fuel, Oxidizer
import numpy as np
import matplotlib.pyplot as plt

methane = Fuel("CH4", mols=1)
o2 = Oxidizer("O2", mols=1)
n2 = Oxidizer("N2", mols=3.76)


# You can specify problems like this!
p1 = Rocket_Problem(pressure=14.7, phi=1, filename="problem1", massf=True)
data = p1.run_cea(methane, o2, n2)

# Or like this
p1 = Rocket_Problem(pressure=14.7, phi=1, materials=[methane, o2, n2], filename="problem1", massf=True)
data = p1.run_cea()
print("Data from stoichiometric combustion of methane/air")

# And if you want to change things
p1 = Rocket_Problem(pressure=14.7, phi=1, materials=[methane, o2, n2], filename="problem1", massf=True)
p1.set_pressure(1000)
p1.set_o_f(17)
p1.set_filename("problem1.1")
p1.set_materials([methane, o2])
data = p1.run_cea()
print("Data from methane/air combustion for o_f=17")
print(data)

# Neat data collector
# let's iterate over pressures!
p1 = Rocket_Problem(phi=1, materials=[methane, o2, n2], filename="problem1", massf=True)

# Collector to get cstar, ideal isp, chamber specific heat ratio, exit density, and throat density as well as fraction of H2O in exhaust
collector = Data_Collector("cstar", "isp", "c_gamma", "c_t", "rho", "t_rho", exit_keys=["H2O"]) # one of chamber_keys or exit_keys, not both
collector1 = Data_Collector(chamber_keys=["H2O"]) # or the other
pressures = np.linspace(14, 10000)
for pressure in pressures:
  print("Running Pressure:", pressure)
  p1.set_pressure(pressure)
  data = p1.run_cea()
  #print(data)
  #input()
  collector.add_data(data)
  collector1.add_data(data)
  
# Access with dot notation
plt.plot(pressures, collector.cstar)
plt.title("Characterstic Velocity")
plt.ylabel("C* (m/s)")
plt.xlabel("Pressure (psi)")
plt.show()

# Or access with dictionary notation
plt.plot(pressures, collector["c_t"])
plt.title("Chamber Temperature")
plt.ylabel("Temperature (K)")
plt.xlabel("Pressure (psi)")
plt.show()

# Plot exit vs chamber products
plt.plot(pressures, collector.H2O, label="Chamber")
plt.plot(pressures, collector1.H2O, label="Exit")
plt.title("Chamber vs Exit H2O mass fraction")
plt.show()
