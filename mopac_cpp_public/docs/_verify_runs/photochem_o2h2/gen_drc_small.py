import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
# Small DRC feasibility: O2 + 1 excited(stretched) H2 aimed at O + 2 ground H2
# O2 along z (O1 at 0,0,0 ; O2 at 0,0,1.21); excited H2 stretched along z,
# near H at z=2.0 heading toward O2; two ground H2 off to the sides.
content = ("PM7 DRC KINE=50\nH2/O2 mix, 1 excited (stretched) H2, DRC dynamics\n\n"
           "O  0.000  0.000  0.000\n"
           "O  0.000  0.000  1.210\n"
           "H  0.000  0.000  2.000\n"
           "H  0.000  0.000  3.500\n"
           "H  2.000  0.000  0.000\n"
           "H  2.000  0.000  0.740\n"
           "H -2.000  0.000  0.000\n"
           "H -2.000  0.000  0.740\n")
with open(os.path.join(D, "drc_small.mop"), "w", newline="\n") as f:
    f.write(content)
print("wrote drc_small.mop")
