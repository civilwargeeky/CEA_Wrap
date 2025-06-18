"""
Tests for CEA_Wrap Problem classes.

These tests focus on ensuring consistent outputs from Problem.run() methods
when given identical inputs. Tests use nested-dictionary based comparisons
to be compatible with future Output class changes.
"""

import unittest
import os
from CEA_Wrap import Fuel, Oxidizer, DetonationProblem, HPProblem, TPProblem, RocketProblem

from common import _convert_to_json_serializable, get_filepath

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
                filename=get_filepath("test_det_simple")
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
                filename=get_filepath("test_det_methane")
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
                filename=get_filepath("test_hp_simple")
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
                filename=get_filepath("test_hp_methane")
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
                filename=get_filepath("test_hp_massf")
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
                filename=get_filepath("test_tp_simple")
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
                filename=get_filepath("test_tp_units")
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
                filename=get_filepath("test_rocket_simple")
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
                filename=get_filepath("test_rocket_solid")
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
                filename=get_filepath("test_rocket_frozen")
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
                filename=get_filepath("test_rocket_pip")
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
                filename=get_filepath("test_rocket_fac")
            )
        
        self._run_consistency_test(create_problem, "Rocket with finite area combustor")


class TestProblemParameterVariation(TestProblemConsistency):
    """Test that problems produce different but consistent results when parameters change."""
    
    def test_pressure_variation_consistency(self):
        """Test that changing pressure produces different but consistent results."""
        materials = [self.h2_gas, self.o2_gas]
        
        # Create problems with different pressures
        problem_low = HPProblem(pressure=500, materials=materials, phi=1.0, filename=get_filepath("test_var_low"))
        problem_high = HPProblem(pressure=2000, materials=materials, phi=1.0, filename=get_filepath("test_var_high"))
        
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


class TestProblemAgainstReferences(TestProblemConsistency):
    """Test that Problem.run() produces outputs consistent with reference JSON files."""
    
    def setUp(self):
        """Set up common materials and load reference directory."""
        super().setUp()
        from pathlib import Path
        self.ref_dir = Path("test_reference_outputs")
        
        # Check if reference directory exists
        if not self.ref_dir.exists():
            self.skipTest("Reference outputs directory not found. Run create_dry_run_outputs() first.")
    
    def _load_reference(self, test_name):
        """Load reference JSON data for a test case."""
        import json
        json_file = self.ref_dir / f"{test_name}.json"
        
        if not json_file.exists():
            self.skipTest(f"Reference file {json_file} not found.")
        
        with open(json_file, "r") as f:
            return json.load(f)
    
    def _compare_with_reference(self, actual_output, reference_data, tolerance=1e-10):
        """Compare actual output with reference data using nested dictionary comparison."""
        # Convert actual output to dictionary
        actual_dict = _convert_to_json_serializable(actual_output)
        
        # Compare the dictionaries
        self._compare_nested_dicts(actual_dict, reference_data, tolerance)
    
    def _compare_nested_dicts(self, actual, reference, tolerance=1e-10, path=""):
        """Recursively compare nested dictionaries with tolerance."""
        # Check that both have the same keys
        actual_keys = set(actual.keys()) if isinstance(actual, dict) else set()
        ref_keys = set(reference.keys()) if isinstance(reference, dict) else set()
        
        # if actual_keys != ref_keys:
        #     missing_in_actual = ref_keys - actual_keys
        #     missing_in_ref = actual_keys - ref_keys
        #     msg = f"Key mismatch at {path}: "
        #     if missing_in_actual:
        #         msg += f"missing in actual: {missing_in_actual} "
        #     if missing_in_ref:
        #         msg += f"missing in reference: {missing_in_ref}"
        #     self.fail(msg)
        
        if (missing_in_actual := ref_keys - actual_keys):
            msg = f"Key mismatch at {path}: "
            msg += f"missing in actual: {missing_in_actual} "
            self.fail(msg)
        
        # Compare values for each key
        for key in ref_keys:
            current_path = f"{path}.{key}" if path else key
            actual_val = actual[key]
            ref_val = reference[key]
            
            if isinstance(actual_val, dict) and isinstance(ref_val, dict):
                # Recursively compare nested dictionaries
                self._compare_nested_dicts(actual_val, ref_val, tolerance, current_path)
            elif isinstance(actual_val, (int, float)) and isinstance(ref_val, (int, float)):
                # Compare numeric values with tolerance
                self.assertAlmostEqual(float(actual_val), float(ref_val), delta=tolerance,
                                     msg=f"Numeric values differ at {current_path}: {actual_val} vs {ref_val}")
            else:
                # Compare non-numeric values exactly
                self.assertEqual(actual_val, ref_val,
                               msg=f"Values differ at {current_path}: {actual_val} vs {ref_val}")
    
    def _test_against_reference(self, test_name, problem_factory):
        """Helper method to test a problem against its reference output."""
        # Load reference data
        reference_data = self._load_reference(test_name)
        
        # Create and run the problem
        problem = problem_factory()
        result = problem.run()
        
        # Compare with reference
        self._compare_with_reference(result, reference_data)
    
    # Sample reference tests - these will test against the JSON files
    def test_hp_h2o2_p500psi_reference(self):
        """Test H2/O2 HP problem at 500 psi against reference."""
        def create_problem():
            return HPProblem(
                pressure=500, materials=[self.h2_gas, self.o2_gas], 
                phi=1.0, pressure_units="psi", filename=get_filepath("test_ref_hp_h2o2_p500psi")
            )
        
        self._test_against_reference("hp_h2o2_p500psi", create_problem)
    
    def test_rocket_h2o2_p1000psi_reference(self):
        """Test H2/O2 rocket problem at 1000 psi against reference."""
        def create_problem():
            return RocketProblem(
                pressure=1000, materials=[self.h2_liquid, self.o2_liquid], 
                phi=1.0, sup=10, pressure_units="psi", filename=get_filepath("test_ref_rocket_h2o2_p1000psi")
            )
        
        self._test_against_reference("rocket_h2o2_p1000psi", create_problem)
    
    def test_det_ch4air_p10bar_reference(self):
        """Test CH4/air detonation at 10 bar against reference."""
        def create_problem():
            return DetonationProblem(
                pressure=10, materials=[self.ch4, self.air], 
                phi=1.0, pressure_units="bar", filename=get_filepath("test_ref_det_ch4air_p10bar")
            )
        
        self._test_against_reference("det_ch4air_p10bar", create_problem)
    
    def test_hp_ch4air_phi1_0_reference(self):
        """Test CH4/air HP at phi=1.0 against reference."""
        def create_problem():
            return HPProblem(
                pressure=1000, materials=[self.ch4, self.air], 
                phi=1.0, filename=get_filepath("test_ref_hp_ch4air_phi1_0")
            )
        
        self._test_against_reference("hp_ch4air_phi1.0", create_problem)
    
    def test_rocket_h2o2_sup5_reference(self):
        """Test H2/O2 rocket with sup=5 against reference."""
        def create_problem():
            return RocketProblem(
                pressure=2000, materials=[self.h2_gas, self.o2_gas], 
                phi=1.0, sup=5, filename=get_filepath("test_ref_rocket_h2o2_sup5")
            )
        
        self._test_against_reference("rocket_h2o2_sup5.0", create_problem)


if __name__ == "__main__":
    unittest.main()