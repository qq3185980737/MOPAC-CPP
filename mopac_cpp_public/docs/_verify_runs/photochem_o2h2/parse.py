import re, os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
names = ["h2_excited","h2_excitedE","o2_gap_bir","o2_singlet2","o2_triplet","o2_singlet"]
for n in names:
    p = os.path.join(D, n+".out")
    if not os.path.exists(p):
        print(f"== {n}: NO OUT"); continue
    txt = open(p, encoding="latin-1").read()
    lines = txt.splitlines()
    emp = ""
    for ln in lines:
        if "Empirical Formula" in ln:
            emp = ln.strip(); break
    heat = tot = ""
    for i,ln in enumerate(lines):
        if "FINAL HEAT OF FORMATION" in ln: heat = ln.strip()
        if "TOTAL ENERGY" in ln and "STATE" in ln: tot = ln.strip()
    print(f"===== {n} =====")
    print("  emp:", emp)
    print("  heat:", heat)
    print("  total:", tot)
    # state table
    for i,ln in enumerate(lines):
        if "STATE       ENERGY" in ln and "RELATIVE" not in ln:
            print("  -- state table --")
            for j in range(i+1, min(i+16, len(lines))):
                s = lines[j]
                if s.strip()=="" or "'RELATIVE" in s or "'The" in s or "'This" in s:
                    continue
                if s.strip().startswith(("1+","2+","3+","4+","5+","6+","1-","2-","3-")):
                    print("    ", s.strip())
            break
