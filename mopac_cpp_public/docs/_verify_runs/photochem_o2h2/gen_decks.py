import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
tests = {
 "o2_triplet.mop": ("PM7 UHF MS=3\no2 triplet\n\nO 0.000 0.000 0.000\nO 0.000 0.000 1.208\n"),
 "o2_singlet.mop": ("PM7 MECI C.I.=2 MS=1\no2 singlet 1Delta\n\nO 0.000 0.000 0.000\nO 0.000 0.000 1.208\n"),
 "h2_excited.mop": ("PM7 MECI C.I.=2 MS=1\nh2 excited\n\nH 0.000 0.000 0.000\nH 0.000 0.000 0.740\n"),
}
for name, content in tests.items():
    with open(os.path.join(D, name), "w", newline="\n") as f:
        f.write(content)
    print("wrote", name)
