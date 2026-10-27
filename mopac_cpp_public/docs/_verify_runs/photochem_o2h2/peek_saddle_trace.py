import re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
lines = open(D+"/ref_ohh2sad.out", encoding="latin-1", errors="replace").read().splitlines()
pts = []
i = 0
n = len(lines)
while i < n:
    s = lines[i].strip()
    m = re.match(r"FOR POINT\s+(\d+)\s+(FIRST|SECOND)\s+STRUCTURE", s)
    if m:
        no = int(m.group(1)); which = m.group(2)
        # find DISTANCE A - B, HEAT, GRADIENT NORM in next few lines
        dist = heat = grad = None
        for j in range(i+1, min(i+6, n)):
            t = lines[j].strip()
            md = re.search(r"DISTANCE A - B\s+([-\d.]+)", t)
            mh = re.search(r"HEAT\s+([-\d.]+)", t)
            mg = re.search(r"GRADIENT NORM\s+([-\d.]+)", t)
            if md: dist = float(md.group(1))
            if mh: heat = float(mh.group(1))
            if mg: grad = float(mg.group(1))
        pts.append((no, which, dist, heat, grad))
        i += 1
    else:
        i += 1
print("total FOR POINT:", len(pts))
if pts:
    mx = max(pts, key=lambda p: p[3])
    print("MAX heat point:", mx)
    first = [p for p in pts if p[1]=="FIRST"]
    second = [p for p in pts if p[1]=="SECOND"]
    print("FIRST:", len(first), "minHeat", min(p[3] for p in first), "maxHeat", max(p[3] for p in first), "first2", [(p[0],p[3]) for p in first[:2]], "last", [(p[0],p[3]) for p in first[-1:]])
    print("SECOND:", len(second), "minHeat", min(p[3] for p in second), "maxHeat", max(p[3] for p in second), "first2", [(p[0],p[3]) for p in second[:2]], "last", [(p[0],p[3]) for p in second[-1:]])
