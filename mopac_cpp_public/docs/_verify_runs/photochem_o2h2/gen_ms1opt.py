import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
content = ("PM7 UHF MS=1\nO2 true triplet S=1 geometry optimization\n\n"
           "O 0.000 0.000 0.000\nO 0.000 0.000 1.208\n")
with open(os.path.join(D, "o2_triplet_opt.mop"), "w", newline="\n") as f:
    f.write(content)
print("wrote o2_triplet_opt.mop")
