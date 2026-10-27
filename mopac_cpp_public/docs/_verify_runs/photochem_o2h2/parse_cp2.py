import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
for name in ["cp2"]:
    txt = open(os.path.join(D, name+".out"), encoding="latin-1").read()
    lines = txt.splitlines()
    emp = next((x.strip() for x in lines if "Empirical Formula" in x), "")
    alpha = next((x.strip() for x in lines if "ALPHA ELECTRONS" in x), "")
    heat = next((x.strip() for x in lines if "FINAL HEAT OF FORMATION" in x), "")
    tot = next((x.strip() for x in lines if "TOTAL ENERGY" in x), "")
    print("==", name)
    print(" emp:", emp)
    print(" alpha:", alpha)
    print(" heat:", heat)
    print(" total:", tot)
