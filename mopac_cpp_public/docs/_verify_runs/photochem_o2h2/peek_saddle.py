import re, os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
n = "ref_ohh2sad"
p = D + "/" + n + ".out"
print("exists:", os.path.exists(p), "len:", os.path.getsize(p) if os.path.exists(p) else None)
if not os.path.exists(p): raise SystemExit
txt = open(p, encoding="latin-1", errors="replace").read()
lines = txt.splitlines()
print("total lines:", len(lines))
print("--- tail ---")
for ln in lines[-18:]:
    print("  ", ln.rstrip())
print("--- markers ---")
for ln in lines:
    s = ln.strip()
    if re.search(r"DISTANCE A - B|HEAT|GRADIENT NORM|FINAL HEAT|TOTAL ENERGY|JOB ENDED|BEST ESTIMATE|C1=|TRANSITION|FIRST|SECOND|SAME SIDE|ERROR|TERMINAT", s, re.I):
        print("  ", s[:110])
