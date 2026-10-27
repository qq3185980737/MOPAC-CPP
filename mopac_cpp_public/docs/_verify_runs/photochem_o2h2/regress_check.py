import re, os
D = "D:/test_rep/mopac_cpp_public/docs/_verify_runs/photochem_o2h2"
def getval(n, pat):
    p = D + "/" + n + ".out"
    if not os.path.exists(p): return "MISSING"
    txt = open(p, encoding="latin-1", errors="replace").read()
    m = re.search(pat, txt)
    return m.group(1) if m else "??"
cases = [
    ("cpp_o2opt",  r"FINAL HEAT OF FORMATION\s*=\s*(-?[\d.]+)", "O2 MS=1 EF HEAT"),
    ("cpp_o2opt",  r"TOTAL ENERGY\s*=\s*(-?[\d.]+)",            "O2 MS=1 EF TOTAL"),
    ("cpp_o2ms1",  r"FINAL HEAT OF FORMATION\s*=\s*(-?[\d.]+)", "O2 MS=1 1SCF HEAT"),
    ("cpp_o2ms1",  r"TOTAL ENERGY\s*=\s*(-?[\d.]+)",            "O2 MS=1 1SCF TOTAL"),
    ("cpp_h2E",    r"FINAL HEAT OF FORMATION\s*=\s*(-?[\d.]+)", "H2 EXCITED HEAT"),
    ("cpp_h2E",    r"TOTAL ENERGY\s*=\s*(-?[\d.]+)",            "H2 EXCITED TOTAL"),
]
ref = {"O2 MS=1 EF HEAT":"-9.17145","O2 MS=1 EF TOTAL":"-587.42567",
       "O2 MS=1 1SCF HEAT":"-4.59425","O2 MS=1 1SCF TOTAL":"-587.22720",
       "H2 EXCITED HEAT":"228.95372","H2 EXCITED TOTAL":"-16.73052"}
for n, pat, label in cases:
    got = getval(n, pat)
    print(f"{label:22s} C++={got:>16s}  ref={ref[label]:>16s}  {'OK' if got==ref[label] else '<<< DIVERGE'}")
