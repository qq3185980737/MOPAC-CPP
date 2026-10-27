import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
# SADDLE: two geometries, same 4 atoms in same order:
#   atom1=O, atom2=H_oh (stays on O), atom3=H_a (transfers H2->O), atom4=H_b (leaves)
# Reactant: OH + H2 (separated along +z)
# Product : H2O + H (H_a bonded to O; H_b free)
content = ("PM7 SADDLE UHF\nOH + H2 -> H2O + H saddle (SADDLE two-geometry)\n\n"
           "O  0.000  0.000  0.000\n"
           "H  0.000  0.000  0.970\n"
           "H  0.000  0.000  3.000\n"
           "H  0.000  0.000  3.740\n"
           "\n"
           "O  0.000  0.000  0.000\n"
           "H  0.000 -0.757  0.586\n"
           "H  0.000  0.757  0.586\n"
           "H  0.000  0.000  3.500\n")
with open(os.path.join(D, "ohh2_saddle.mop"), "w", newline="\n") as f:
    f.write(content)
print("wrote ohh2_saddle.mop")
