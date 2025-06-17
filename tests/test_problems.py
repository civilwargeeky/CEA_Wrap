"""
Tests for CEA_Wrap Problem classes.

These tests focus on ensuring consistent outputs from Problem.run() methods
when given identical inputs. Tests use nested-dictionary based comparisons
to be compatible with future Output class changes.
"""

import unittest
from CEA_Wrap import Fuel, Oxidizer, DetonationProblem, HPProblem, TPProblem, RocketProblem


class TestProblemConsistency(unittest.TestCase):
    """Test that Problem.run() produces consistent outputs for identical inputs."""
    
    def setUp(self):
        """Set up common materials for testing."""
        # Common fuels
        self.h2_liquid = Fuel("H2(L)", temp=20)
        self.h2_gas = Fuel("H2")
        self.ch4 = Fuel("CH4")
        self.aluminum = Fuel("AL(cr)")
        self.htpb = Fuel("HTPB")
        
        # Common oxidizers
        self.o2_liquid = Oxidizer("O2(L)", temp=90)
        self.o2_gas = Oxidizer("O2")
        self.air = Oxidizer("Air")
        self.ap = Oxidizer("NH4CLO4(I)")
        self.h2o2 = Oxidizer("H2O2(L)")
    
    def _assert_output_consistency(self, output1, output2, tolerance=1e-10):
        """
        Assert that two Output objects are consistent within tolerance.
        Uses nested dictionary comparison to be compatible with future Output changes.
        """
        # Convert outputs to dictionaries for comparison
        dict1 = dict(output1)
        dict2 = dict(output2)
        
        # Check that both have the same keys
        self.assertEqual(set(dict1.keys()), set(dict2.keys()),
                        "Output objects have different keys")
        
        # Check values for each key
        for key in dict1.keys():
            val1, val2 = dict1[key], dict2[key]
            
            # Handle nested Output objects (like prod_c, prod_e, etc.)
            if hasattr(val1, '__dict__') and hasattr(val2, '__dict__'):
                self._assert_output_consistency(val1, val2, tolerance)
            elif isinstance(val1, (int, float)) and isinstance(val2, (int, float)):
                self.assertAlmostEqual(val1, val2, delta=tolerance,
                                     msg=f"Values for key '{key}' differ: {val1} vs {val2}")
            else:
                self.assertEqual(val1, val2, 
                               msg=f"Non-numeric values for key '{key}' differ: {val1} vs {val2}")

    def _run_consistency_test(self, problem_factory, description=""):
        """
        Helper method to test that a problem produces consistent results
        when run multiple times with identical inputs.
        """
        # Create two identical problems and run them
        problem1 = problem_factory()
        problem2 = problem_factory() 
        
        result1 = problem1.run()
        result2 = problem2.run()
        
        # Assert they produce identical results
        self._assert_output_consistency(result1, result2)
        
        # Also test that running the same problem twice gives consistent results
        result3 = problem1.run()
        self._assert_output_consistency(result1, result3)


class TestDetonationProblem(TestProblemConsistency):
    """Tests for DetonationProblem consistency."""
    
    def test_simple_detonation_consistency(self):
        """Test simple H2/O2 detonation problem consistency."""
        def create_problem():
            return DetonationProblem(
                pressure=1000,
                materials=[self.h2_gas, self.o2_gas],
                phi=1.0,
                filename="test_det_simple"
            )
        
        self._run_consistency_test(create_problem, "Simple H2/O2 detonation")
    
    def test_methane_air_detonation_consistency(self):
        """Test methane/air detonation problem consistency."""
        def create_problem():
            return DetonationProblem(
                pressure=500,
                materials=[self.ch4, self.air],
                phi=0.8,
                pressure_units="psi",
                filename="test_det_methane"
            )
        
        self._run_consistency_test(create_problem, "Methane/air detonation")


class TestHPProblem(TestProblemConsistency):
    """Tests for HPProblem (constant enthalpy and pressure) consistency."""
    
    def test_simple_hp_consistency(self):
        """Test simple H2/O2 HP problem consistency."""
        def create_problem():
            return HPProblem(
                pressure=1000,
                materials=[self.h2_gas, self.o2_gas],
                phi=1.0,
                filename="test_hp_simple"
            )
        
        self._run_consistency_test(create_problem, "Simple H2/O2 HP")
    
    def test_methane_air_hp_consistency(self):
        """Test methane/air HP problem consistency."""
        def create_problem():
            return HPProblem(
                pressure=2000,
                materials=[self.ch4, self.air],
                o_f=15.0,
                pressure_units="psi",
                filename="test_hp_methane"
            )
        
        self._run_consistency_test(create_problem, "Methane/air HP")
    
    def test_hp_with_massf_consistency(self):
        """Test HP problem with mass fractions output consistency."""
        def create_problem():
            return HPProblem(
                pressure=1500,
                materials=[self.h2_gas, self.o2_gas],
                phi=1.2,
                massf=True,
                filename="test_hp_massf"
            )
        
        self._run_consistency_test(create_problem, "HP with mass fractions")


class TestTPProblem(TestProblemConsistency):
    """Tests for TPProblem (constant temperature and pressure) consistency."""
    
    def test_simple_tp_consistency(self):
        """Test simple aluminum/air TP problem consistency."""
        def create_problem():
            return TPProblem(
                pressure=1,
                temperature=2350,
                materials=[self.aluminum, self.air],
                phi=1.0,
                pressure_units="atm",
                filename="test_tp_simple"
            )
        
        self._run_consistency_test(create_problem, "Simple Al/air TP")
    
    def test_tp_different_units_consistency(self):
        """Test TP problem with different unit systems consistency."""
        def create_problem():
            return TPProblem(
                pressure=100,
                temperature=2000,
                materials=[self.aluminum, self.air],
                phi=1.0,
                pressure_units="bar",
                temperature_units="r",
                filename="test_tp_units"
            )
        
        self._run_consistency_test(create_problem, "TP with different units")


class TestRocketProblem(TestProblemConsistency):
    """Tests for RocketProblem consistency."""
    
    def test_simple_rocket_consistency(self):
        """Test simple H2/O2 rocket problem consistency."""
        def create_problem():
            return RocketProblem(
                pressure=2000,
                materials=[self.h2_liquid, self.o2_liquid],
                phi=1.0,
                sup=5,
                filename="test_rocket_simple"
            )
        
        self._run_consistency_test(create_problem, "Simple H2/O2 rocket")
    
    def test_solid_rocket_consistency(self):
        """Test solid propellant rocket problem consistency."""
        def create_problem():
            # Set up solid propellant materials with weight percentages
            al = Fuel("AL(cr)", wt_percent=10)
            htpb = Fuel("HTPB", wt_percent=20)  
            ap = Oxidizer("NH4CLO4(I)", wt_percent=70)
            
            problem = RocketProblem(
                pressure=1000,
                materials=[al, htpb, ap],
                sup=15,
                filename="test_rocket_solid"
            )
            problem.set_absolute_o_f()  # Calculate o/f from weight percentages
            return problem
        
        self._run_consistency_test(create_problem, "Solid rocket")
    
    def test_frozen_rocket_consistency(self):
        """Test frozen flow rocket problem consistency."""
        def create_problem():
            return RocketProblem(
                pressure=2000,
                materials=[self.h2_liquid, self.o2_liquid],
                phi=1.0,
                sup=5,
                analysis_type="frozen",
                nfz=2,  # Freeze at throat
                filename="test_rocket_frozen"
            )
        
        self._run_consistency_test(create_problem, "Frozen rocket")
    
    def test_rocket_with_pip_consistency(self):
        """Test rocket problem with pressure ratio consistency."""
        def create_problem():
            return RocketProblem(
                pressure=1000,
                materials=[self.ch4, self.o2_gas],
                o_f=3.5,
                pip=68,  # Pressure ratio instead of area ratio
                filename="test_rocket_pip"
            )
        
        self._run_consistency_test(create_problem, "Rocket with pressure ratio")
    
    def test_rocket_with_fac_consistency(self):
        """Test rocket with finite area combustor consistency."""
        def create_problem():
            return RocketProblem(
                pressure=1500,
                materials=[self.h2_gas, self.o2_gas],
                phi=1.0,
                sup=10,
                fac_ac=2.0,
                filename="test_rocket_fac"
            )
        
        self._run_consistency_test(create_problem, "Rocket with finite area combustor")


class TestProblemParameterVariation(TestProblemConsistency):
    """Test that problems produce different but consistent results when parameters change."""
    
    def test_pressure_variation_consistency(self):
        """Test that changing pressure produces different but consistent results."""
        materials = [self.h2_gas, self.o2_gas]
        
        # Create problems with different pressures
        problem_low = HPProblem(pressure=500, materials=materials, phi=1.0, filename="test_var_low")
        problem_high = HPProblem(pressure=2000, materials=materials, phi=1.0, filename="test_var_high")
        
        result_low = problem_low.run()
        result_high = problem_high.run()
        
        # Results should be different (temperature should be different)
        self.assertNotAlmostEqual(result_low.t, result_high.t, delta=1.0,
                                msg="Temperature should differ significantly with pressure change")
        
        # But running each again should give consistent results
        result_low_2 = problem_low.run()
        result_high_2 = problem_high.run()
        
        self._assert_output_consistency(result_low, result_low_2)
        self._assert_output_consistency(result_high, result_high_2)


def create_dry_run_outputs():
    """
    Utility function to generate reference outputs for all test cases.
    Run this once to create expected outputs, then use those for comparison.
    """
    import json
    import pickle
    from pathlib import Path
    
    # Create reference directory
    ref_dir = Path("test_reference_outputs")
    ref_dir.mkdir(exist_ok=True)
    
    # Define all test cases
    test_cases = {
        "det_simple": lambda: DetonationProblem(
            pressure=1000, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=1.0, filename="ref_det_simple"
        ).run(),
        
        "hp_simple": lambda: HPProblem(
            pressure=1000, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=1.0, filename="ref_hp_simple"
        ).run(),
        
        "tp_simple": lambda: TPProblem(
            pressure=1, temperature=2350, materials=[Fuel("AL(cr)"), Oxidizer("Air")], 
            phi=1.0, pressure_units="atm", filename="ref_tp_simple"
        ).run(),
        
        "rocket_simple": lambda: RocketProblem(
            pressure=2000, materials=[Fuel("H2(L)", temp=20), Oxidizer("O2(L)", temp=90)], 
            phi=1.0, sup=5, filename="ref_rocket_simple"
        ).run(),
    }
    
    # Generate and save reference outputs
    references = {}
    for name, problem_func in test_cases.items():
        try:
            result = problem_func()
            # Convert to dictionary for serialization
            references[name] = dict(result)
            print(f"Generated reference for {name}")
        except Exception as e:
            print(f"Failed to generate reference for {name}: {e}")
    
    # Save references
    with open(ref_dir / "references.pkl", "wb") as f:
        pickle.dump(references, f)
    
    print(f"Saved {len(references)} reference outputs to {ref_dir}")


if __name__ == "__main__":
    # Uncomment the line below to generate reference outputs
    # create_dry_run_outputs()
    
    unittest.main()