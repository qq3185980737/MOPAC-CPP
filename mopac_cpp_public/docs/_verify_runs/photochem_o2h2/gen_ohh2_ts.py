import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
# OH + H2 -> H2O + H, near-linear saddle guess
# O(0,0,0); Ha(0,0,1.35) forming O-H; Hb(0,0,2.18) leaving H (breaking Ha-Hb);
# H_oh at 0.97 A from O, ~104 deg from O-Ha axis
content = ("PM7 TS UHF PRECISE\nOH + H2 -> H2O + H saddle point (linear guess)\n\n"
           "O  0.000  0.000  0.000\n"
           "H  0.000  0.000  1.350\n"
           "H  0.000  0.000  2.180\n"
           "H  0.941  0.000 -0.235\n")
with open(os.path.join(D, "ohh2_ts.mop"), "w", newline="\n") as f:
    f.write(content)
print("wrote ohh2_ts.mop")
