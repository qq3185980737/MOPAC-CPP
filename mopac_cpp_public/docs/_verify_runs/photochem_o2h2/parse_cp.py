import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
name = "cp1"
txt = open(os.path.join(D, name+".out"), encoding="latin-1").read()
lines = txt.splitlines()
emp = next((x.strip() for x in lines if "Empirical Formula" in x), "")
heat = next((x.strip() for x in lines if "FINAL HEAT OF FORMATION" in x), "")
tot = next((x.strip() for x in lines if "TOTAL ENERGY" in x and "STATE" in x), "")
print("==", name)
print(" emp:", emp)
print(" heat:", heat)
print(" total:", tot)
for i, x in enumerate(lines):
    if "STATE       ENERGY" in x and "RELATIVE" not in x:
        print(" --state table--")
        for j in range(i+1, min(i+16, len(lines))):
            s = lines[j].strip()
            if s and not (s.startswith("'") or s.startswith("RELATIVE") or s.startswith("The") or s.startswith("This")):
                print("   ", s)
        break
