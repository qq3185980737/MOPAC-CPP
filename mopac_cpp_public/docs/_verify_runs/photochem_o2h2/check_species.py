import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
for s in ['H','O','OH','H2','O2','H2O']:
    for tag in ['ref','cpp']:
        p = os.path.join(D, f"{tag}_sp{s}.out")
        print(f"{tag}_sp{s}: exists={os.path.exists(p)} size={os.path.getsize(p) if os.path.exists(p) else None}")
# show a few failing outputs' tails
for s in ['H','O','OH','H2','O2']:
    p = os.path.join(D, f"ref_sp{s}.out")
    if os.path.exists(p):
        print(f"\n--- ref_sp{s}.out tail ---")
        t = open(p, encoding="latin-1", errors="replace").read().splitlines()
        for ln in t[-8:]:
            print("  ", ln.rstrip())
