import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
tests = {
 "t2_diff.mop": ("PM7 UHF MS=2\nHO test\nH 0.0 0.0 0.0\nO 0.0 0.0 1.0\n"),
 "t2_same.mop": ("PM7 UHF MS=3\ndiatomic same element\nO 0.0 0.0 0.0\nO 0.0 0.0 1.0\n"),
 "t3_same.mop": ("PM7 UHF MS=1\nthree O atoms\nO 0.0 0.0 0.0\nO 0.0 0.0 1.0\nO 2.0 0.0 1.0\n"),
}
for name, content in tests.items():
    with open(os.path.join(D, name), "w", newline="\n") as f:
        f.write(content)
    print("wrote", name)
