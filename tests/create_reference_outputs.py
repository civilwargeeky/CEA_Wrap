from CEA_Wrap import Fuel, Oxidizer, DetonationProblem, HPProblem, TPProblem, RocketProblem

from common import get_filepath, _convert_to_json_serializable

def create_dry_run_outputs():
    """
    Utility function to generate extensive reference outputs for all test cases.
    Run this once to create expected outputs, then use those for comparison.
    This creates a comprehensive matrix of test cases covering various parameter combinations.
    
    Each test case is saved as a separate JSON file in the test_reference_outputs directory.
    """
    import json
    from pathlib import Path
    import os
    
    # Create reference directory
    ref_dir = Path("test_reference_outputs")
    ref_dir.mkdir(exist_ok=True)
    
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
            phi=1.0, pressure_units="psi", filename=get_filepath(f"ref_hp_p{p}psi")
        ).run()
        
        test_cases[f"rocket_h2o2_p{p}psi"] = lambda p=p: RocketProblem(
            pressure=p, materials=[Fuel("H2(L)", temp=20), Oxidizer("O2(L)", temp=90)], 
            phi=1.0, sup=10, pressure_units="psi", filename=get_filepath(f"ref_rocket_p{p}psi")
        ).run()
    
    for p in pressures_bar:
        test_cases[f"det_ch4air_p{p}bar"] = lambda p=p: DetonationProblem(
            pressure=p, materials=[Fuel("CH4"), Oxidizer("Air")], 
            phi=1.0, pressure_units="bar", filename=get_filepath(f"ref_det_p{p}bar")
        ).run()
    
    # ================================================================
    # O/F AND PHI RATIO VARIATIONS
    # ================================================================
    phi_ratios = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0]
    of_ratios = [2.0, 4.0, 8.0, 15.0, 20.0]
    
    for phi in phi_ratios:
        test_cases[f"hp_ch4air_phi{phi}"] = lambda phi=phi: HPProblem(
            pressure=1000, materials=[Fuel("CH4"), Oxidizer("Air")], 
            phi=phi, filename=get_filepath(f"ref_hp_phi{phi}")
        ).run()
        
        test_cases[f"rocket_h2o2_phi{phi}"] = lambda phi=phi: RocketProblem(
            pressure=2000, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=phi, sup=15, filename=get_filepath(f"ref_rocket_phi{phi}")
        ).run()
    
    for of in of_ratios:
        test_cases[f"hp_ch4air_of{of}"] = lambda of=of: HPProblem(
            pressure=1500, materials=[Fuel("CH4"), Oxidizer("Air")], 
            o_f=of, filename=get_filepath(f"ref_hp_of{of}")
        ).run()
        
        test_cases[f"det_h2o2_of{of}"] = lambda of=of: DetonationProblem(
            pressure=1000, materials=[Fuel("H2"), Oxidizer("O2")], 
            o_f=of, filename=get_filepath(f"ref_det_of{of}")
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
                phi=1.0, filename=get_filepath(f"ref_hp_ft{ft}_ot{ot}")
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
            phi=1.0, sup=20, filename=get_filepath(f"ref_rocket_cryo_{ft}_{ot}")
        ).run()
    
    # ================================================================
    # ROCKET AREA RATIO VARIATIONS (SUP/SUB)
    # ================================================================
    sup_ratios = [1.5, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    sub_ratios = [10, 3, 2, 1.5, 1.1]
    
    for sup in sup_ratios:
        test_cases[f"rocket_h2o2_sup{sup}"] = lambda sup=sup: RocketProblem(
            pressure=2000, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=1.0, sup=sup, filename=get_filepath(f"ref_rocket_sup{sup}")
        ).run()
    
    for sub in sub_ratios:
        test_cases[f"rocket_ch4air_sub{sub}"] = lambda sub=sub: RocketProblem(
            pressure=1000, materials=[Fuel("CH4"), Oxidizer("Air")], 
            phi=1.0, sub=sub, filename=get_filepath(f"ref_rocket_sub{sub}")
        ).run()
    
    # Pressure ratio (pip) variations
    pip_ratios = [2, 5, 10, 20, 50, 100]
    for pip in pip_ratios:
        test_cases[f"rocket_h2o2_pip{pip}"] = lambda pip=pip: RocketProblem(
            pressure=2000, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=1.0, pip=pip, filename=get_filepath(f"ref_rocket_pip{pip}")
        ).run()
    
    # ================================================================
    # FINITE AREA COMBUSTOR (FAC) VARIATIONS
    # ================================================================
    fac_ac_ratios = [1.2, 1.5, 2.0, 3.0, 5.0]
    fac_ma_ratios = [0.1, 0.3, 0.5, 0.7]
    
    for fac in fac_ac_ratios:
        test_cases[f"rocket_h2o2_fac_ac{fac}"] = lambda fac=fac: RocketProblem(
            pressure=1500, materials=[Fuel("H2"), Oxidizer("O2")], 
            phi=1.0, sup=10, fac_ac=fac, filename=get_filepath(f"ref_rocket_fac_ac{fac}")
        ).run()
    
    for fac in fac_ma_ratios:
        test_cases[f"rocket_ch4o2_fac_ma{fac}"] = lambda fac=fac: RocketProblem(
            pressure=1000, materials=[Fuel("CH4"), Oxidizer("O2")], 
            phi=1.0, sup=15, fac_ma=fac, filename=get_filepath(f"ref_rocket_fac_ma{fac}")
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
            filename=get_filepath(f"ref_rocket_frozen_nfz{nfz}")
        ).run()
    
    for custom_nfz in custom_nfz_ratios:
        test_cases[f"rocket_ch4air_frozen_custom{custom_nfz}"] = lambda custom_nfz=custom_nfz: RocketProblem(
            pressure=1500, materials=[Fuel("CH4"), Oxidizer("Air")], 
            phi=1.0, sup=25, analysis_type="frozen", custom_nfz=custom_nfz, 
            filename=get_filepath(f"ref_rocket_frozen_custom{custom_nfz}")
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
            filename=get_filepath(f"ref_tp_t{temp}k")
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
            filename=get_filepath(f"ref_tp_units_{temp_unit}_{press_unit}")
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
                sup=15, filename=get_filepath(f"ref_solid_{name}")
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
        phi=1.8, sup=50, filename=get_filepath("ref_rocket_complex_high_p_rich")
    ).run()
    
    # Low pressure, fuel-lean, subsonic
    test_cases["rocket_complex_low_p_lean"] = lambda: RocketProblem(
        pressure=200, materials=[Fuel("CH4", temp=400), Oxidizer("Air", temp=600)], 
        phi=0.6, sub=3.0, filename=get_filepath("ref_rocket_complex_low_p_lean")
    ).run()
    
    # Monopropellant test
    test_cases["rocket_monoprop_h2o2"] = lambda: RocketProblem(
        pressure=1000, materials=[Fuel("H2O2(L)")], 
        o_f=0, sup=10, filename=get_filepath("ref_rocket_monoprop_h2o2")
    ).run()
    
    # Mass fraction outputs
    test_cases["hp_massf_complex"] = lambda: HPProblem(
        pressure=1500, materials=[Fuel("CH4", temp=600), Oxidizer("Air", temp=800)], 
        phi=1.3, massf=True, filename=get_filepath("ref_hp_massf_complex")
    ).run()
    
    # ================================================================
    # GENERATE AND SAVE ALL REFERENCE OUTPUTS  
    # ================================================================
    
    successful_cases = []
    failed_cases = []
    
    total_cases = len(test_cases)
    print(f"Generating {total_cases} reference test cases...")
    
    for i, (name, problem_func) in enumerate(test_cases.items(), 1):
        try:
            print(f"[{i:3d}/{total_cases}] Generating {name}...")
            result = problem_func()
            
            # Convert to JSON-serializable format
            json_data = _convert_to_json_serializable(result)
            
            # Save individual JSON file
            json_file = ref_dir / f"{name}.json"
            with open(json_file, "w") as f:
                json.dump(json_data, f, indent=2)
            
            successful_cases.append(name)
            
        except Exception as e:
            print(f"[{i:3d}/{total_cases}] FAILED {name}: {e}")
            failed_cases.append((name, str(e)))
    
    # Save a summary report
    with open(ref_dir / "_generation_report.txt", "w") as f:
        f.write(f"Reference Generation Report\n")
        try:
            from importlib.metadata import version
            f.write(f"CEA Wrap Version: {version('CEA_Wrap')}\n")
        except:
            f.write(f"CEA Wrap Version: Unknown\n")
        f.write(f"==========================\n\n")
        f.write(f"Total test cases: {total_cases}\n")
        f.write(f"Successful: {len(successful_cases)}\n")
        f.write(f"Failed: {len(failed_cases)}\n\n")
        
        if successful_cases:
            f.write("Successful Cases:\n")
            f.write("-----------------\n")
            for name in successful_cases:
                f.write(f"{name}.json\n")
            f.write("\n")
        
        if failed_cases:
            f.write("Failed Cases:\n")
            f.write("-------------\n")
            for name, error in failed_cases:
                f.write(f"{name}: {error}\n")
    
    print(f"\nGeneration complete!")
    print(f"Saved {len(successful_cases)} successful reference outputs as individual JSON files to {ref_dir}")
    if failed_cases:
        print(f"{len(failed_cases)} cases failed - see _generation_report.txt for details")

if __name__ == "__main__":
    create_dry_run_outputs()