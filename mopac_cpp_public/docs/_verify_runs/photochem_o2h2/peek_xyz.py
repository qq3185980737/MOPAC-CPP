import os, math, re
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
p = D+"/ref_drcs.xyz"
raw = open(p, encoding="latin-1", errors="replace").read()
lines = raw.splitlines()
# DRC xyz: frames of natoms+2 lines? or plain xyz natoms/comment/atoms.
# Let's detect structure
print("total lines:", len(lines))
print("first 12:")
for l in lines[:12]:
    print("  |", l.rstrip())
print("last 8:")
for l in lines[-8:]:
    print("  |", l.rstrip())
