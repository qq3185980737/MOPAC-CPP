import os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
def mk(name, kw, comment, atoms):
    body = "\n".join(f"{lab} {x:.3f} {y:.3f} {z:.3f}" for lab,x,y,z in atoms)
    return (f"{kw}\n{comment}\n\n{body}\n")
tests = {
 "h2_excited.mop":     mk("h2_excited.mop","PM7 MECI C.I.=2 MS=1","h2 excited via CI", [("H",0,0,0),("H",0,0,0.74)]),
 "h2_excitedE.mop":    mk("h2_excitedE.mop","PM7 MECI C.I.=2 EXCITED","h2 excited state E", [("H",0,0,0),("H",0,0,0.74)]),
 "o2_gap_bir.mop":     mk("o2_gap_bir.mop","PM7 UHF MECI C.I.=2 BIRAD","o2 biradical CI", [("O",0,0,0),("O",0,0,1.208)]),
 "o2_singlet2.mop":    mk("o2_singlet2.mop","PM7 MECI C.I.=2 OPEN=2","o2 open singlet", [("O",0,0,0),("O",0,0,1.208)]),
}
for name, content in tests.items():
    with open(os.path.join(D, name), "w", newline="\n") as f:
        f.write(content)
    print("wrote", name)
