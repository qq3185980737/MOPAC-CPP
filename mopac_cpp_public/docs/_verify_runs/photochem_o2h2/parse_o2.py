import os, re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
for n in ["o2_triplet","cpp_o2t"]:
    p = os.path.join(D, n+".out")
    if not os.path.exists(p) or os.path.getsize(p)==0:
        print(f"{n}: NO/EMPTY OUT"); continue
    txt = open(p, encoding="latin-1", errors="replace").read()
    lines = txt.splitlines()
    print(f"===== {n} =====")
    pats = ["Empirical Formula","ALPHA ELECTRONS","FINAL HEAT OF FORMATION","TOTAL ENERGY","STATE:","MOLECULAR WEIGHT","DGEMM","ERROR","Crashed","JOB ENDED","COMPUTATION TIME","1SCF","GRADIENT"]
    for ln in lines:
        s = ln.strip()
        if any(k in s for k in pats):
            print("  ", s)
    # tail
    print("  -- tail --")
    for s in lines[-8:]:
        print("  ", s.strip())
