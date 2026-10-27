import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
def w(name, body):
    with open(os.path.join(D, name+".mop"), "w", newline="\n") as f:
        f.write(body)
# species optimized heats: PM7; UHF for open shell, MS=1 for triplet O/O2
w("sp_H",  "PM7 UHF\nH atom (doublet)\n\nH  0.0 0.0 0.0\n")
w("sp_O",  "PM7 UHF MS=1\nO atom (triplet 3P)\n\nO  0.0 0.0 0.0\n")
w("sp_OH", "PM7 UHF\nOH radical (doublet)\n\nO  0.0 0.0 0.0\nH  0.0 0.0 0.970\n")
w("sp_H2", "PM7\nH2 (singlet)\n\nH  0.0 0.0 0.0\nH  0.0 0.0 0.740\n")
w("sp_O2", "PM7 UHF MS=1\nO2 (triplet)\n\nO  0.0 0.0 0.0\nO  0.0 0.0 1.208\n")
w("sp_H2O","PM7\nH2O (singlet)\n\nO  0.0 0.0 0.0\nH  0.0 0.757 0.586\nH  0.0 -0.757 0.586\n")
print("wrote 6 species decks")
