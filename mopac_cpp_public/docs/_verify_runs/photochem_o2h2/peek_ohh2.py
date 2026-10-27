import re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
txt = open(D+"/ref_ohh2ts.out", encoding="latin-1", errors="replace").read()
lines = txt.splitlines()
print("total lines:", len(lines))
print("--- tail ---")
for ln in lines[-16:]:
    print("  ", ln.rstrip())
print("--- markers ---")
for ln in lines:
    s = ln.strip()
    if re.search(r"FINAL HEAT|TOTAL ENERGY|GRADIENT|CYCLE|JOB|IMAGINAR|FREQUENC|MODE|EIGENVAL|1ST |saddle|TRANSITION|ORDER|NEGATIVE", s, re.I):
        print("  ", s[:110])
