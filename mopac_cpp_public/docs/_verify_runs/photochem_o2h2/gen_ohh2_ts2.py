import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
# OH + H2 -> H2O + H, literature PM7 saddle guess (bent O-Ha-Hb ~170deg)
# O(0,0,0); Ha(0,0,1.345) forming O-H; Hb=(0.143,0,2.158) leaving (Ha-Hb=0.825);
# H_oh=(0.941,0,-0.166) O-H=0.955, Hoh-O-Ha~100deg
content = ("PM7 TS UHF PRECISE GEO-OK\nOH + H2 -> H2O + H saddle (bent guess)\n\n"
           "O  0.000  0.000  0.000\n"
           "H  0.000  0.000  1.345\n"
           "H  0.143  0.000  2.158\n"
           "H  0.941  0.000 -0.166\n")
with open(os.path.join(D, "ohh2_ts2.mop"), "w", newline="\n") as f:
    f.write(content)
print("wrote ohh2_ts2.mop")
