import os, subprocess, re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
REF = "D:/test_rep/mopac_fortran_ref/mopac_ref.exe"
os.chdir(D)
elem = ["O","O","H","H","H","H"]
base = [
 [3.000,4.000,4.000],
 [3.000,4.000,5.210],
 [5.500,4.000,4.000],
 [5.500,4.000,4.740],
 [3.000,6.500,4.000],
 [3.000,6.500,4.740],
]
nat = len(base)

def deck(geo):
    lines = ["PM7 C.I.=2 EXCITED 1SCF", "grad check", ""]
    for i in range(nat):
        lines.append(f"{elem[i]}  {geo[i][0]:.4f}  {geo[i][1]:.4f}  {geo[i][2]:.4f}")
    return "\n".join(lines) + "\n"

def run_hof(name, geo):
    with open(name+".mop","w",newline="\n") as f: f.write(deck(geo))
    subprocess.run([REF, name+".mop"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t=open(name+".out",encoding="latin-1",errors="replace").read()
    m=re.search(r"FINAL HEAT OF FORMATION\s*=\s*([-\d.]+)",t)
    return float(m.group(1)) if m else None

# Analytic gradient parsed from EF cycle-1 GRADIENTS block
analyt = [
  1.500,   1.500, -109.581,    # O1
  1.303,   1.303,  109.176,    # O2
 -1.381,   0.022,   27.589,    # H1
 -1.466,   0.022,  -27.386,    # H2
  0.022,  -1.381,   27.589,    # H3
  0.022,  -1.466,  -27.386,    # H4
]
delta = 0.01
E0 = run_hof("gn_base", base)
print(f"E0(base) = {E0} kcal/mol")
numer = []
for a in range(nat):
    for c in range(3):
        gp = [r[:] for r in base]; gm = [r[:] for r in base]
        gp[a][c]+=delta; gm[a][c]-=delta
        Ep = run_hof(f"gn_p_{a}_{c}", gp)
        Em = run_hof(f"gn_m_{a}_{c}", gm)
        g = (Ep-Em)/(2*delta)
        numer.append(g)
        print(f"  atom{a+1}({elem[a]}) coord{c+1}: Ep={Ep:.5f} Em={Em:.5f}  numeric_g={g:+.3f}  analytic_g={analyt[3*a+c]:+.3f}  diff={g-analyt[3*a+c]:+.3f}")
print("\nSUM |diff| =", sum(abs(numer[i]-analyt[i]) for i in range(3*nat)))
