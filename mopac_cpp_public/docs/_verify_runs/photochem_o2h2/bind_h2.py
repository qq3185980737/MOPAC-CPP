import os, re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
def total(tag, s):
    t = open(os.path.join(D, f"{tag}_sp{s}.out"), encoding="latin-1", errors="replace").read()
    m = re.search(r"TOTAL ENERGY\s*=\s*(-?[\d.]+)", t)
    return float(m.group(1)) if m else None
for tag in ['ref','cpp']:
    eH = total(tag,"H"); eH2 = total(tag,"H2")
    bind = eH2 - 2*eH   # negative = bound
    bindkcal = bind * 23.0605
    dhfH = 52.10200
    print(f"{tag}: E(H)={eH}  E(H2)={eH2}  binding={bind:.4f} eV = {bindkcal:.3f} kcal/mol  -> D(H2) expected ~-4.5eV")
