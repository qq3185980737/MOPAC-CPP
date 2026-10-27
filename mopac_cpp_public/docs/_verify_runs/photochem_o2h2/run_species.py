import os, re, subprocess
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
REF = r"D:\test_rep\mopac_fortran_ref\mopac_ref.exe"
CPP = r"D:\test_rep\mopac_cpp_public\x64\Release\mopac.exe"
species = ['H','O','OH','H2','O2','H2O']
def run(exe, name):
    r = subprocess.run([exe, name], cwd=D, capture_output=True)
    return r.returncode
def heat(name):
    p = os.path.join(D, name+".out")
    if not os.path.exists(p): return "NONE"
    t = open(p, encoding="latin-1", errors="replace").read()
    m = re.search(r"FINAL HEAT OF FORMATION\s*=\s*(-?[\d.]+)", t)
    return m.group(1) if m else "??"
for s in species:
    for exe, tag in [(REF,"ref"),(CPP,"cpp")]:
        n = f"{tag}_sp{s}"
        shutil = None
        import shutil as sh
        sh.copyfile(os.path.join(D,f"sp_{s}.mop"), os.path.join(D,n))
        for ext in (".out",".arc"):
            pp = os.path.join(D,n+ext)
            if os.path.exists(pp): os.remove(pp)
        ec = run(exe, n)
    print(f"{s:4s} ref_heat={heat('ref_sp'+s):>14s}  cpp_heat={heat('cpp_sp'+s):>14s}")
