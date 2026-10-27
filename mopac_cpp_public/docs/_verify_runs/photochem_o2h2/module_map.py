import os, re, glob

f90dir = "D:/mopac/mopac-archive/2016_MOPAC2016/MOPAC_source_code"
cppdir = "D:/test_rep/mopac_cpp_public"

# Fortran module source files (routines) — .F90 not interface (_I) or .mod
f90 = [os.path.basename(p) for p in glob.glob(os.path.join(f90dir, "*.F90"))]
f90 = [f for f in f90 if not f.endswith("_I.F90") and not f.endswith(".F90.bak")]

cpp = [os.path.basename(p) for p in glob.glob(os.path.join(cppdir, "*.cpp"))]
cph = [os.path.basename(p) for p in glob.glob(os.path.join(cppdir, "*.h"))]

def stem(name):
    return re.sub(r"(_C)?\.(F90|f90|cpp|h)$", "", name).lower()

# map f90 stem -> has cpp or h in cpp tree (by stem match)
missing = []
present = []
for f in sorted(f90):
    s = stem(f)
    hit_cpp = any(stem(c)==s for c in cpp) or any(stem(h)==s for h in cph)
    if hit_cpp:
        present.append((f, s))
    else:
        missing.append((f, s))

print("TOTAL Fortran modules:", len(f90))
print("== PRESENT in C++ ({}) ==".format(len(present)))
print(", ".join(s for _,s in present))
print("== MISSING in C++ ({}) ==".format(len(missing)))
for f,s in missing:
    print("  ", f)
