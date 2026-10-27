import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
kcal = 23.0605
# PM7 species heats (ref == cpp, bit-identical except H2O ~2e-5)
H, O, OH, H2, O2, H2O = 52.10200, 59.55900, 7.25699, -32.01058, -9.17145, -57.79968
# physical (experimental, kcal/mol): H2~0, O2~0, H2O~-57.8, OH~9.3, H~52.1, O~59.6
phys = {"H2":0.0, "O2":0.0, "H2O":-57.8, "OH":9.3}
print("="*78)
print("PHOTOLYSIS-INITIATED H2+O2 CHAIN: PM7 heats (ref==cpp) vs experimental")
print("="*78)
print(f"{'species':6s} {'PM7 dHf':>10s} {'exp dHf':>9s} {'note'}")
print(f"{'H':6s} {H:10.3f} {52.10:9.1f}")
print(f"{'O':6s} {O:10.3f} {59.56:9.1f}")
print(f"{'OH':6s} {OH:10.3f} {9.3:9.1f}")
print(f"{'H2':6s} {H2:10.3f} {0.0:9.1f}  <<< PM7 overbinds H2 by ~32 kcal (both programs agree)")
print(f"{'O2':6s} {O2:10.3f} {0.0:9.1f}")
print(f"{'H2O':6s} {H2O:10.3f} {-57.8:9.1f}")
print()
steps = [
    ("INIT (photolysis):  H2 + hv -> 2H",        2*H - H2),
    ("INIT (photolysis):  O2 + hv -> 2O",         2*O - O2),
    ("chain:  H  + O2 -> OH + O",                 (OH+O) - (H+O2)),
    ("chain:  O  + H2 -> OH + H",                 (OH+H) - (O+H2)),
    ("chain:  OH + H2 -> H2O + H",                (H2O+H) - (OH+H2)),
    ("chain:  H  + OH -> H2O",                    H2O - (H+OH)),
    ("OVERALL: 2H2 + O2 -> 2H2O",                 2*H2O - (2*H2+O2)),
]
print("STEP DELTA-H (PM7, kcal/mol):")
print("-"*78)
for name, dh in steps:
    ev = dh/kcal
    print(f"  {name:38s} dH = {dh:+9.3f} kcal = {ev:+7.3f} eV")
print("-"*78)
print("NOTE: every step consuming H2 (INIT H2, O+H2, OH+H2, OVERALL) is shifted")
print("by +~32 kcal (endothermic-er) because PM7 overbinds H2 by ~32 kcal/mol.")
