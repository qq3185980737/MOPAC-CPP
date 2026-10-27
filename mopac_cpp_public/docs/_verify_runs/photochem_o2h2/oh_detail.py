import os, math
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
            pr = lines[i].split()
            atoms.append((pr[0], float(pr[1]), float(pr[2]), float(pr[3])))
            i += 1
        # extract heat from comment
        import re
        m = re.search(r"HEAT OF FORMATION =\s*([-\d.]+)", comment)
        heat = float(m.group(1)) if m else None
        frames.append((comment, atoms, heat))
    else:
        i += 1
nat = len(frames[0][1])
# atoms: 0,1 = O ; 2,3 = excited H2 (H3,H4) ; 4,5 = ground H2a ; 6,7 = ground H2b
print("frames:", len(frames), "natoms:", nat)
# focus on the O2 (atoms 0,1) and excited-H2 H's (2,3): track O-H distances first 300 fs
print("\n--- first 300 fs: O(any)-H(excited H2) distances, heat ---")
print(" fs    O0-H3   O0-H4   O1-H3   O1-H4   heat")
for fi in range(0, 300, 4):
    if fi >= len(frames): break
    c,a,h = frames[fi]
    fs = fi*0.1
    d = {}
    for ai in (2,3):
        for oj in (0,1):
            d[(oj,ai)] = math.dist(a[oj][1:], a[ai][1:])
    print(f" {fs:6.1f} {d[(0,2)]:6.3f} {d[(0,3)]:6.3f} {d[(1,2)]:6.3f} {d[(1,3)]:6.3f}  {h if h is not None else 0:8.3f}")
# min O-H for excited H2 over whole run + frame
best = None
for fi,(c,a,h) in enumerate(frames):
    dmin = min(math.dist(a[o][1:], a[hh][1:]) for o in (0,1) for hh in (2,3))
    if best is None or dmin < best[1]: best = (fi, dmin)
print("\nmin O-H (O2 vs excited-H2 H's): frame", best[0], "d=", round(best[1],3), "fs=", round(best[0]*0.1,1))
# residence below 1.6 A
below = [fi for fi,(c,a,h) in enumerate(frames) if min(math.dist(a[o][1:],a[hh][1:]) for o in (0,1) for hh in (2,3)) < 1.6]
print("frames with O(excitedH)-H < 1.6 A:", len(below), "range", (below[0],below[-1]) if below else None)
