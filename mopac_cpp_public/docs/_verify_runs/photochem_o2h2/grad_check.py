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

def deck(geo, kw):
    lines = [kw, "grad check", ""]
    for i in range(nat):
        lines.append(f"{elem[i]}  {geo[i][0]:.4f}  {geo[i][1]:.4f}  {geo[i][2]:.4f}")
    return "\n".join(lines) + "\n"

def run(name, geo, kw):
    with open(name+".mop", "w", newline="\n") as f:
        f.write(deck(geo, kw))
    subprocess.run([REF, name+".mop"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t = open(name+".out", encoding="latin-1", errors="replace").read()
    return t

def hof(t):
    m = re.search(r"FINAL HEAT OF FORMATION\s*=\s*([-\d.]+)", t)
    return float(m.group(1)) if m else None

# 1) Analytic gradient (Fortran DERIV prints GRADIENTS)
t = run("gr_a", base, "PM7 C.I.=2 EXCITED 1SCF DERIV")
print("analytic E0 =", hof(t))
# find GRADIENTS block
i = t.find("GRADIENTS")
seg = t[i:i+400] if i>=0 else ""
nums = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", seg)
analyt = [float(x) for x in nums[:3*nat]] if len(nums)>=3*nat else None
if analyt:
    print("analytic grad (kcal/mol/A):")
    for i2 in range(nat):
        print(f"  atom{i2+1}({elem[i2]}): " + " ".join(f"{analyt[3*i2+c]:+.4f}" for c in range(3)))
else:
    print("!! analytic gradient not parsed; search markers:")
    print("  GRADIENTS found at", i)
    print(seg[:250])
