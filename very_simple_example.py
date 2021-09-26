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

# And if you want to change things
p1 = Rocket_Problem(pressure=14.7, phi=1, materials=[methane, o2, n2], filename="problem1", massf=True)
p1.set_pressure(1000)
p1.set_o_f(17)
p1.set_filename("problem1.1")
p1.set_materials([methane, o2])
data = p1.run_cea()

# Neat data collector
# let's iterate over pressures!
p1 = Rocket_Problem(phi=1, materials=[methane, o2, n2], filename="problem1", massf=True)

collector = Data_Collector("cstar", "isp", "gamma", "c_t", "rho", chamber_keys=["H2O"]) # one 
collector = Data_Collector("cstar", "isp", "gamma", "c_t", "rho", exit_keys=["H2O"]) # or the other
pressures = np.linspace(14, 10000)
for pressure in pressures:
  print("Running Pressure:", pressure)
  p1.set_pressure(pressure)
  data = p1.run_cea()
  #print(data)
  #input()
  collector.add_data(data)
  

plt.plot(pressures, collector.cstar)
plt.show()

plt.plot(pressures, collector["c_t"])
plt.show()

plt.plot(pressures, collector.H2O)
plt.show()