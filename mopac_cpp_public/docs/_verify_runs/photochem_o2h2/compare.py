import sys, os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"

def extract(path):
    txt = open(path, encoding="latin-1", errors="replace").read()
    lines = txt.splitlines()
    out = {}
    for i, ln in enumerate(lines):
        s = ln.strip()
        if "FINAL HEAT OF FORMATION" in s: out["heat"] = s
        if "TOTAL ENERGY" in s and "STATE" in s: out["total"] = s
        if "Empirical Formula" in s and "emp" not in out: out["emp"] = s
        if "STATE       ENERGY" in s and "RELATIVE" not in s:
            st = []
            for j in range(i+1, min(i+16, len(lines))):
                t = lines[j].strip()
                if t == "" or t.startswith("'"): continue
                if t[:2] in ("1+","2+","3+","4+","5+","1-","2-","3-") or (t and t[0].isdigit() and "+" in t[:4]):
                    st.append(t)
            out["states"] = st
        if "COMPONENT OF SPIN" in s and "spin" not in out: out["spin"] = s
        if "NO OF CONFIGURATIONS" in s and "nconf" not in out: out["nconf"] = s
    return out

a = extract(os.path.join(D, "h2_excitedE.out"))   # reference
b = extract(os.path.join(D, "cpp_h2E.out"))       # c++ port
print("FIELD          REFERENCE                                    C++ PORT")
for k in ("emp","heat","total","spin","nconf"):
    av = a.get(k, "-"); bv = b.get(k, "-")
    print(f"{k:6s}  {av:50s}  {bv}")
print("\n--- REF state table ---")
for s in a.get("states", []): print("  ", s)
print("--- C++ state table ---")
for s in b.get("states", []): print("  ", s)
