import os, math, re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
p = D+"/ref_drcs.xyz"
lines = open(p, encoding="latin-1", errors="replace").read().splitlines()
frames = []
i = 0
while i < len(lines):
    if lines[i].strip().isdigit():
        n = int(lines[i].strip()); i += 1
        comment = lines[i]; i += 1
        atoms = []
        for k in range(n):
            parts = lines[i].split()
            atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            i += 1
        frames.append((comment, atoms))
    else:
        i += 1
print("frames:", len(frames))
def ohmin(atoms):
    best = 1e9; pair=None
    O=[a for a in atoms if a[0]=='O']; H=[a for a in atoms if a[0]=='H']
    for a in O:
        for b in H:
            d = math.dist(a[1:], b[1:])
            if d < best: best, pair = d, (a[0],b[0])
    return best
mins = []
for ci,(comment, atoms) in enumerate(frames):
    if ci % 20 != 0: continue
    m = ohmin(atoms)
    mins.append((ci, m))
# find global min
best=None
for ci,(comment,atoms) in enumerate(frames):
    m=ohmin(atoms)
    if best is None or m<best[1]: best=(ci,m,atoms)
print("min O-H over sampled frames:", [(ci,round(m,3)) for ci,m in mins])
print("global min O-H: frame", best[0], "d=", round(best[1],3))
if best[1] < 1.5:
    print("BEST FRAME atoms:")
    for a in best[2]: print("  ", a)
print("O-H dist at final frame:", round(ohmin(frames[-1][1]),3))
print("O-O dists:")
for ci,(comment,atoms) in enumerate([frames[0], frames[len(frames)//2], frames[-1]]):
    O=[a for a in atoms if a[0]=='O']
    dd=[round(math.dist(O[0][1:],O[1][1:]),3)] if len(O)==2 else []
    print("  frame", ci, "O-O", dd)
