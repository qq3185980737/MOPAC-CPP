import os, re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
for n in ["ref_o2opt","cpp_o2opt"]:
    txt = open(os.path.join(D, n+".out"), encoding="latin-1", errors="replace").read()
    lines = txt.splitlines()
    print(f"===== {n} =====")
    for ln in lines:
        s = ln.strip()
        if re.search(r"FINAL HEAT|TOTAL ENERGY|ELECTRONIC ENERGY|ALPHA ELECTRONS|BETA ELECTRONS|Empirical|FINAL GEOMETRY|O         |MOLECULAR WEIGHT|CARTESIAN|GRADIENT", s):
            print("  ", s[:120])
