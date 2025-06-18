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
    Utility function to generate extensive reference outputs for all test cases.
    Run this once to create expected outputs, then use those for comparison.
    This creates a comprehensive matrix of test cases covering various parameter combinations.
    """
    import json
    import pickle
    from pathlib import Path
    import os
    
    # Create reference directory
    ref_dir = Path("test_reference_outputs")
    ref_dir.mkdir(exist_ok=True)

    int_dir = "test_reference_interstitial"
    if not os.path.isdir(int_dir):
        os.makedirs(int_dir)
    
    # Define comprehensive test cases
    test_cases = {}
    
    # ================================================================
    # PRESSURE VARIATIONS
    # ================================================================
    pressures_psi = [500, 1000, 2000, 3000]
    pressures_bar = [10, 50, 100, 200]
    pressures_atm = [1, 10, 50, 100]
    
    for p in pressures_psi:
        test_cases[f"hp_h2o2_p{p}psi"] = lambda p=p: HPProblem(
            pressure=p, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=1.0, pressure_units="psi", filename=int_dir+os.sep+f"ref_hp_p{p}psi"
        ).run()
        
        test_cases[f"rocket_h2o2_p{p}psi"] = lambda p=p: RocketProblem(
            pressure=p, materials=[Fuel("H2(L)", temp=20), Oxidizer("O2(L)", temp=90)], 
            phi=1.0, sup=10, pressure_units="psi", filename=int_dir+os.sep+f"ref_rocket_p{p}psi"
        ).run()
    
    for p in pressures_bar:
        test_cases[f"det_ch4air_p{p}bar"] = lambda p=p: DetonationProblem(
            pressure=p, materials=[Fuel("CH4"), Oxidizer("Air")], 
            phi=1.0, pressure_units="bar", filename=int_dir+os.sep+f"ref_det_p{p}bar"
        ).run()
    
    # ================================================================
    # O/F AND PHI RATIO VARIATIONS
    # ================================================================
    phi_ratios = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0]
    of_ratios = [2.0, 4.0, 8.0, 15.0, 20.0]
    
    for phi in phi_ratios:
        test_cases[f"hp_ch4air_phi{phi}"] = lambda phi=phi: HPProblem(
            pressure=1000, materials=[Fuel("CH4"), Oxidizer("Air")], 
            phi=phi, filename=int_dir+os.sep+f"ref_hp_phi{phi}"
        ).run()
        
        test_cases[f"rocket_h2o2_phi{phi}"] = lambda phi=phi: RocketProblem(
            pressure=2000, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=phi, sup=15, filename=int_dir+os.sep+f"ref_rocket_phi{phi}"
        ).run()
    
    for of in of_ratios:
        test_cases[f"hp_ch4air_of{of}"] = lambda of=of: HPProblem(
            pressure=1500, materials=[Fuel("CH4"), Oxidizer("Air")], 
            o_f=of, filename=int_dir+os.sep+f"ref_hp_of{of}"
        ).run()
        
        test_cases[f"det_h2o2_of{of}"] = lambda of=of: DetonationProblem(
            pressure=1000, materials=[Fuel("H2"), Oxidizer("O2")], 
            o_f=of, filename=int_dir+os.sep+f"ref_det_of{of}"
        ).run()
    
    # ================================================================
    # TEMPERATURE VARIATIONS
    # ================================================================
    fuel_temps = [200, 298, 500, 800, 1000]
    ox_temps = [200, 298, 400, 600, 800]
    
    for ft in fuel_temps:
        for ot in ox_temps[::2]:  # Sample every other ox temp to reduce combinations
            test_cases[f"hp_h2o2_ft{ft}_ot{ot}"] = lambda ft=ft, ot=ot: HPProblem(
                pressure=1000, 
                materials=[Fuel("H2", temp=ft), Oxidizer("O2", temp=ot)], 
                phi=1.0, filename=int_dir+os.sep+f"ref_hp_ft{ft}_ot{ot}"
            ).run()
    
    # Cryogenic temperatures
    cryo_cases = [
        (20, 90, "H2(L)", "O2(L)"),
        (110, 90, "CH4(L)", "O2(L)"),
        (298, 298, "N2H4", "N2O4(L)")
    ]
    
    for ft, ot, fname, oname in cryo_cases:
        test_cases[f"rocket_cryo_{fname.replace('(', '').replace(')', '').lower()}_{oname.replace('(', '').replace(')', '').lower()}"] = lambda ft=ft, ot=ot, fname=fname, oname=oname: RocketProblem(
            pressure=1500, 
            materials=[Fuel(fname, temp=ft), Oxidizer(oname, temp=ot)], 
            phi=1.0, sup=20, filename=int_dir+os.sep+f"ref_rocket_cryo_{ft}_{ot}"
        ).run()
    
    # ================================================================
    # ROCKET AREA RATIO VARIATIONS (SUP/SUB)
    # ================================================================
    sup_ratios = [1.5, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    sub_ratios = [10, 3, 2, 1.5, 1.1]
    
    for sup in sup_ratios:
        test_cases[f"rocket_h2o2_sup{sup}"] = lambda sup=sup: RocketProblem(
            pressure=2000, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=1.0, sup=sup, filename=int_dir+os.sep+f"ref_rocket_sup{sup}"
        ).run()
    
    for sub in sub_ratios:
        test_cases[f"rocket_ch4air_sub{sub}"] = lambda sub=sub: RocketProblem(
            pressure=1000, materials=[Fuel("CH4"), Oxidizer("Air")], 
            phi=1.0, sub=sub, filename=int_dir+os.sep+f"ref_rocket_sub{sub}"
        ).run()
    
    # Pressure ratio (pip) variations
    pip_ratios = [2, 5, 10, 20, 50, 100]
    for pip in pip_ratios:
        test_cases[f"rocket_h2o2_pip{pip}"] = lambda pip=pip: RocketProblem(
            pressure=2000, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=1.0, pip=pip, filename=int_dir+os.sep+f"ref_rocket_pip{pip}"
        ).run()
    
    # ================================================================
    # FINITE AREA COMBUSTOR (FAC) VARIATIONS
    # ================================================================
    fac_ac_ratios = [1.2, 1.5, 2.0, 3.0, 5.0]
    fac_ma_ratios = [0.1, 0.3, 0.5, 0.7]
    
    for fac in fac_ac_ratios:
        test_cases[f"rocket_h2o2_fac_ac{fac}"] = lambda fac=fac: RocketProblem(
            pressure=1500, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=1.0, sup=10, fac_ac=fac, filename=int_dir+os.sep+f"ref_rocket_fac_ac{fac}"
        ).run()
    
    for fac in fac_ma_ratios:
        test_cases[f"rocket_ch4o2_fac_ma{fac}"] = lambda fac=fac: RocketProblem(
            pressure=1000, materials=[Fuel("CH4"), Oxidizer("O2")], 
            phi=1.0, sup=15, fac_ma=fac, filename=int_dir+os.sep+f"ref_rocket_fac_ma{fac}"
        ).run()
    
    # ================================================================
    # FROZEN FLOW ANALYSIS VARIATIONS
    # ================================================================
    nfz_positions = [1, 2, 3]  # Chamber, throat, exit
    custom_nfz_ratios = [1.2, 2.0, 5.0, 10.0]
    
    for nfz in nfz_positions:
        test_cases[f"rocket_h2o2_frozen_nfz{nfz}"] = lambda nfz=nfz: RocketProblem(
            pressure=2000, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=1.0, sup=20, analysis_type="frozen", nfz=nfz, 
            filename=int_dir+os.sep+f"ref_rocket_frozen_nfz{nfz}"
        ).run()
    
    for custom_nfz in custom_nfz_ratios:
        test_cases[f"rocket_ch4air_frozen_custom{custom_nfz}"] = lambda custom_nfz=custom_nfz: RocketProblem(
            pressure=1500, materials=[Fuel("CH4"), Oxidizer("Air")], 
            phi=1.0, sup=25, analysis_type="frozen", custom_nfz=custom_nfz, 
            filename=int_dir+os.sep+f"ref_rocket_frozen_custom{custom_nfz}"
        ).run()
    
    # ================================================================
    # TP PROBLEM VARIATIONS
    # ================================================================
    temperatures_k = [1500, 2000, 2500, 3000, 3500]
    temp_units = ["k", "r", "c", "f"]
    
    for temp in temperatures_k:
        test_cases[f"tp_alcr_air_t{temp}k"] = lambda temp=temp: TPProblem(
            pressure=1, temperature=temp, 
            materials=[Fuel("AL(cr)"), Oxidizer("Air")], 
            phi=1.0, pressure_units="atm", temperature_units="k",
            filename=int_dir+os.sep+f"ref_tp_t{temp}k"
        ).run()
    
    # Different unit combinations for TP
    unit_combinations = [
        (2000, "r", 100, "bar"),
        (1500, "c", 50, "psi"),
        (2500, "f", 5, "atm")
    ]
    
    for temp, temp_unit, press, press_unit in unit_combinations:
        test_cases[f"tp_units_{temp_unit}_{press_unit}"] = lambda temp=temp, temp_unit=temp_unit, press=press, press_unit=press_unit: TPProblem(
            pressure=press, temperature=temp, 
            materials=[Fuel("AL(cr)"), Oxidizer("Air")], 
            phi=1.0, pressure_units=press_unit, temperature_units=temp_unit,
            filename=int_dir+os.sep+f"ref_tp_units_{temp_unit}_{press_unit}"
        ).run()
    
    # ================================================================
    # SOLID PROPELLANT VARIATIONS
    # ================================================================
    solid_compositions = [
        (10, 20, 70, "basic_solid"),  # Al/HTPB/AP
        (15, 15, 70, "high_al"),      # Higher aluminum
        (5, 25, 70, "high_htpb"),     # Higher HTPB
        (12, 18, 70, "balanced"),     # Balanced composition
    ]
    
    for al_pct, htpb_pct, ap_pct, name in solid_compositions:
        def create_solid_problem(al_pct=al_pct, htpb_pct=htpb_pct, ap_pct=ap_pct, name=name):
            al = Fuel("AL(cr)", wt_percent=al_pct)
            htpb = Fuel("HTPB", wt_percent=htpb_pct)
            ap = Oxidizer("NH4CLO4(I)", wt_percent=ap_pct)
            problem = RocketProblem(
                pressure=1000, materials=[al, htpb, ap], 
                sup=15, filename=int_dir+os.sep+f"ref_solid_{name}"
            )
            problem.set_absolute_o_f()
            return problem.run()
        
        test_cases[f"rocket_solid_{name}"] = create_solid_problem
    
    # ================================================================
    # COMBINATION TEST CASES (COMPLEX SCENARIOS)
    # ================================================================
    
    # High pressure, fuel-rich, supersonic
    test_cases["rocket_complex_high_p_rich"] = lambda: RocketProblem(
        pressure=5000, materials=[Fuel("H2", temp=500), Oxidizer("O2", temp=300)], 
        phi=1.8, sup=50, filename=int_dir+os.sep+"ref_rocket_complex_high_p_rich"
    ).run()
    
    # Low pressure, fuel-lean, subsonic
    test_cases["rocket_complex_low_p_lean"] = lambda: RocketProblem(
        pressure=200, materials=[Fuel("CH4", temp=400), Oxidizer("Air", temp=600)], 
        phi=0.6, sub=3.0, filename=int_dir+os.sep+"ref_rocket_complex_low_p_lean"
    ).run()
    
    # Monopropellant test
    test_cases["rocket_monoprop_h2o2"] = lambda: RocketProblem(
        pressure=1000, materials=[Fuel("H2O2(L)")], 
        o_f=0, sup=10, filename=int_dir+os.sep+"ref_rocket_monoprop_h2o2"
    ).run()
    
    # Mass fraction outputs
    test_cases["hp_massf_complex"] = lambda: HPProblem(
        pressure=1500, materials=[Fuel("CH4", temp=600), Oxidizer("Air", temp=800)], 
        phi=1.3, massf=True, filename=int_dir+os.sep+"ref_hp_massf_complex"
    ).run()
    
    # ================================================================
    # GENERATE AND SAVE ALL REFERENCE OUTPUTS  
    # ================================================================
    
    references = {}
    failed_cases = []
    
    total_cases = len(test_cases)
    print(f"Generating {total_cases} reference test cases...")
    
    for i, (name, problem_func) in enumerate(test_cases.items(), 1):
        try:
            print(f"[{i:3d}/{total_cases}] Generating {name}...")
            result = problem_func()
            # Convert to dictionary for serialization
            references[name] = dict(result)
        except Exception as e:
            print(f"[{i:3d}/{total_cases}] FAILED {name}: {e}")
            failed_cases.append((name, str(e)))
    
    # Save references
    with open(ref_dir / "references.pkl", "wb") as f:
        pickle.dump(references, f)
    
    # Save a summary report
    with open(ref_dir / "generation_report.txt", "w") as f:
        f.write(f"Reference Generation Report\n")
        f.write(f"==========================\n\n")
        f.write(f"Total test cases: {total_cases}\n")
        f.write(f"Successful: {len(references)}\n")
        f.write(f"Failed: {len(failed_cases)}\n\n")
        
        if failed_cases:
            f.write("Failed Cases:\n")
            f.write("-------------\n")
            for name, error in failed_cases:
                f.write(f"{name}: {error}\n")
    
    print(f"\nGeneration complete!")
    print(f"Saved {len(references)} successful reference outputs to {ref_dir}")
    if failed_cases:
        print(f"{len(failed_cases)} cases failed - see generation_report.txt for details")


if __name__ == "__main__":
    # Uncomment the line below to generate reference outputs
    create_dry_run_outputs()
    
    # unittest.main()