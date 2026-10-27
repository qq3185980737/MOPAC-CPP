import re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
txt = open(D+"/o2_triplet.out", encoding="latin-1", errors="replace").read()
lines = txt.splitlines()
print("== REF MS=3 out total lines:", len(lines))
print("--- tail ---")
for ln in lines[-12:]:
    print("  ", ln.rstrip())
print("--- markers ---")
for ln in lines:
    s = ln.strip()
    if re.search(r"CYCLE|FINAL HEAT|TOTAL ENERGY|dsyevd|CODE|JOB ENDED|GRADIENT", s, re.I):
        print("  ", s[:110])
