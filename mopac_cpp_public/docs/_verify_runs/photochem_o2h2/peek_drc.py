import re, os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
p = D+"/ref_drcs.out"
t = open(p, encoding="latin-1", errors="replace").read()
lines = t.splitlines()
print("total lines:", len(lines))
print("--- DRC markers ---")
for ln in lines:
    s = ln.strip()
    if re.search(r"DRC|VELOCITY|EXCESS KINETIC|TIME|STEP|DYNAMIC|trajector|INTEGRAT|KINETIC ENERGY|POTENTIAL", s, re.I):
        print("  ", s[:100])
        # cap
