import os, re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
for n in ["ref_o2t1","cpp_o2t1"]:
    p = os.path.join(D, n+".out")
    txt = open(p, encoding="latin-1", errors="replace").read()
    lines = txt.splitlines()
    print(f"===== {n} =====")
    for ln in lines:
        s = ln.strip()
        if re.search(r"HEAT|TOTAL ENERGY|ELECTRONIC ENERGY|CORE-CORE|IONIZATION|ALPHA|BETA|Empirical", s):
            print("  ", s)
    print("  tail:", [x.strip() for x in lines[-6:]])
