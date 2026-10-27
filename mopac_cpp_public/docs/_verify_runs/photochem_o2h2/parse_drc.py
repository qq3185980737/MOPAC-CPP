import re, os, math
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
p = D+"/ref_drcs.out"
t = open(p, encoding="latin-1", errors="replace").read()
lines = t.splitlines()
# 1) FEMTOSECONDS table rows
print("--- FEMTOSECONDS trajectory (first 12 / last 6) ---")
fs = [ln for ln in lines if re.search(r"^\s*\d+\.\d+\s+\d+\s", ln) and "TIME" not in ln and "WALL" not in ln]
# better: lines following the FEMTOSECONDS header
idx = [i for i,l in enumerate(lines) if "FEMTOSECONDS" in l]
if idx:
    start = idx[0]
    rows=[]
    for l in lines[start+1:start+4000]:
        s=l.strip()
        if not s: break
        if re.match(r"^[\d\.]+\s", s): rows.append(s)
    print("rows:", len(rows))
    for r in rows[:12]: print("  ", r[:90])
    for r in rows[-6:]: print("  ", r[:90])
# 2) final geometry: find last block of "N  O/H   x y z"
print("--- last full geometry block (final structure) ---")
# find coordinate lines: pattern "   N   O   x y z" near end
all_geom = re.findall(r"^\s*(\d+)\s+([OH])\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$", t, re.M)
print("total atom-coordinate lines found:", len(all_geom))
if all_geom:
    # group into frames of 8
    natoms=8
    nframes = len(all_geom)//natoms
    last = all_geom[-natoms:]
    print("frames:", nframes)
    print("LAST FRAME:")
    for a in last:
        print("  ", a)
    # compute O-H distances in last frame
    atoms=[]
    for a in last:
        atoms.append((a[1], float(a[2]), float(a[3]), float(a[4])))
    print("O-H distances in final frame (<2.0 A):")
    for i,ai in enumerate(atoms):
        for j,aj in enumerate(atoms):
            if i<j and ((ai[0]=='O' and aj[0]=='H') or (ai[0]=='H' and aj[0]=='O')):
                d=math.dist(ai[1:],aj[1:])
                if d<2.0:
                    print(f"   O(i={i}) H(j={j})  d={d:.3f} A")
    print("O-O distance final:", [round(math.dist(a[1:],b[1:]),3) for a in atoms for b in atoms if a[0]=='O' and b[0]=='O' and atoms.index(a)<atoms.index(b)])
