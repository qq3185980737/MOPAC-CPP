import re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
txt = open(D+"/cpp_o2t3fix.out", encoding="latin-1", errors="replace").read()
lines = txt.splitlines()
print("== C++ MS=3 out total lines:", len(lines))
# tail
print("--- tail ---")
for ln in lines[-14:]:
    print("  ", ln.rstrip())
# search markers
print("--- markers ---")
for ln in lines:
    s = ln.strip()
    if re.search(r"ERROR|HEAT|TOTAL ENERGY|DSYEVD|diagonal|CYCLE|GRADIENT|JOB ENDED|abort|failed|converg", s, re.I):
        print("  ", s[:110])
