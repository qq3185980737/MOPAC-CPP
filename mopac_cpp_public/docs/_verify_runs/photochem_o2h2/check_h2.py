import os, re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
for tag in ['ref','cpp']:
    p = os.path.join(D, f"{tag}_spH2.out")
    t = open(p, encoding="latin-1", errors="replace").read()
    print(f"=== {tag}_spH2 ===")
    for ln in t.splitlines():
        s = ln.strip()
        if re.search(r"Empirical Formula|FINAL HEAT|TOTAL ENERGY|CARTESIAN|^\d+\s+H\s| \d+\s*O\s|HEAT OF", s):
            print("  ", s[:100])
