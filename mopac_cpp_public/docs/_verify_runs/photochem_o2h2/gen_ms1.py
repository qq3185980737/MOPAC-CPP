import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
content = ("PM7 UHF MS=1 1SCF\nO2 true triplet (9a 7b) single point\n\n"
           "O 0.000 0.000 0.000\nO 0.000 0.000 1.208\n")
for name in ["o2_triplet_ms1.mop"]:
    with open(os.path.join(D, name), "w", newline="\n") as f:
        f.write(content)
print("wrote", name)
