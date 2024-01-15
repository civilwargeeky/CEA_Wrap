# Quick Example Testing out TP Problems (specified temperature and pressure)

import CEA_Wrap as cw
prob1 = cw.TPProblem(materials=[cw.Fuel("AL(cr)"), cw.Oxidizer("Air")], phi=1, temperature=2350, pressure=1, pressure_units="atm").run()
print(prob1)

# Let's try with different units. Rankine and Bar
prob2 = cw.TPProblem(materials=[cw.Fuel("AL(cr)"), cw.Oxidizer("Air")], phi=1, temperature=2000, temperature_units="r", pressure=100, pressure_units="bar").run()
print(prob2)