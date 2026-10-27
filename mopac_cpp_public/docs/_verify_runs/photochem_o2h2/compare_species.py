import os, re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
sp = ['H','O','OH','H2','O2','H2O']
print("SPECIES HEAT COMPARISON (reference vs C++)")
for s in sp:
    row = []
    for tag in ['ref','cpp']:
        p = os.path.join(D, f"{tag}_sp{s}.out")
        t = open(p, encoding="latin-1", errors="replace").read()
        m = re.search(r"(FINAL |)HEAT OF FORMATION\s*=\s*(-?[\d.]+)", t)
        row.append(m.group(2) if m else "??")
        # also show the raw line for transparency on the first species
    print(f"{s:4s} ref={row[0]:>14s}  cpp={row[1]:>14s}  {'OK' if row[0]==row[1] else ''}")
print("\nraw heat lines:")
for s in ['H','O','OH','H2','O2']:
    for tag in ['ref','cpp']:
        p = os.path.join(D, f"{tag}_sp{s}.out")
        t = open(p, encoding="latin-1", errors="replace").read()
        for ln in t.splitlines():
            if 'HEAT OF FORMATION' in ln:
                print(f"  {tag}_sp{s}: {ln.strip()}")
