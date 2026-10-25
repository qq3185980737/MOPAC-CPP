#!/usr/bin/env python
"""run_moz_opt_probe.py — probe whether MOZYME geometry optimisation is usable
in the C++ port, with and without COSMO, and with several optimisers.

For each (molecule, keyword) pair the C++ port and the official reference are
run and the following are compared:
  exit code, normal termination, final HOF, final GRADIENT NORM,
  number of optimisation cycles, first-cycle gradient, failure messages.
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
CPP = os.path.join(ROOT, "bin", "mopac_cpp.exe")
REF = os.path.join(ROOT, "bin", "mopac_ref.exe")

WATER = """O 0.000 0.000 0.000
H 0.000 0.757 0.586
H 0.000 -0.757 0.586
"""
METHANOL = """C 0.000 0.000 0.000
O 1.430 0.000 0.000
H 1.800 0.890 0.000
H -0.360 1.030 0.000
H -0.360 -0.515 0.890
H -0.360 -0.515 -0.890
"""
BENZENE = """C 1.3970 0.0000 0.0000
C 0.6985 1.2100 0.0000
C -0.6985 1.2100 0.0000
C -1.3970 0.0000 0.0000
C -0.6985 -1.2100 0.0000
C 0.6985 -1.2100 0.0000
H 2.4810 0.0000 0.0000
H 1.2405 2.1490 0.0000
H -1.2405 2.1490 0.0000
H -2.4810 0.0000 0.0000
H -1.2405 -2.1490 0.0000
H 1.2405 -2.1490 0.0000
"""

C12 = """C 0 0 0
H 0 0 1.0
H 0 0 -1.0
H 0 1.0 0
C 1.54 0 0
H 1.54 0 1.0
H 1.54 0 -1.0
C 3.08 0 0
H 3.08 0 1.0
H 3.08 0 -1.0
C 4.62 0 0
H 4.62 0 1.0
H 4.62 0 -1.0
C 6.16 0 0
H 6.16 0 1.0
H 6.16 0 -1.0
C 7.7 0 0
H 7.7 0 1.0
H 7.7 0 -1.0
C 9.24 0 0
H 9.24 0 1.0
H 9.24 0 -1.0
C 10.78 0 0
H 10.78 0 1.0
H 10.78 0 -1.0
C 12.32 0 0
H 12.32 0 1.0
H 12.32 0 -1.0
C 13.86 0 0
H 13.86 0 1.0
H 13.86 0 -1.0
C 15.4 0 0
H 15.4 0 1.0
H 15.4 0 -1.0
C 16.94 0 0
H 16.94 0 1.0
H 16.94 0 -1.0
H 16.94 1.0 0
"""

MOLS = {"water": WATER, "benzene": BENZENE, "c12": C12}
KEYWORDS = [
    ("default", "PM7 MOZYME"),
    ("BFGS", "PM7 MOZYME BFGS"),
    ("L-BFGS", "PM7 MOZYME L-BFGS"),
    ("EF", "PM7 MOZYME EF"),
    ("COSMO", "PM7 MOZYME EPS=78.4"),
]

RE_HOF = re.compile(r"FINAL HEAT OF FORMATION\s*=\s*(-?\d+\.\d+)")
RE_GN = re.compile(r"GRADIENT NORM\s*=\s*(-?\d+\.\d+)")
RE_CYC = re.compile(r"CYCLE:\s*(\d+).*GRAD\.:\s*(-?\d+\.\d+).*HEAT:\s*(-?\d+\.\d+)")
BAD = ["OPTIMIZATION ABORTED", "EXCESSIVE", "FAILED", "NOT CONVERGED",
       "TOO MANY", "JOB ABORTED", "ERROR"]


def run(exe, tag, kw, geo):
    sub = os.path.join(ROOT, "moz_opt_probe", tag, "cpp" if exe == CPP else "ref")
    os.makedirs(sub, exist_ok=True)
    with open(os.path.join(sub, tag + ".mop"), "w") as f:
        f.write(kw + "\n" + tag + "\n\n" + geo)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    with open(os.path.join(sub, "stdout.txt"), "wb") as fo, \
         open(os.path.join(sub, "stderr.txt"), "wb") as fe:
        try:
            p = subprocess.run([exe, tag + ".mop"], cwd=sub, env=env,
                               stdout=fo, stderr=fe, timeout=1800)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            rc = "TIMEOUT"
    text = ""
    out = os.path.join(sub, tag + ".out")
    if os.path.exists(out):
        with open(out, "r", errors="replace") as f:
            text = f.read()
    hof = RE_HOF.findall(text)
    gn = RE_GN.findall(text)
    cyc = RE_CYC.findall(text)
    bad = [b for b in BAD if b in text.upper()]
    return {
        "exit": rc,
        "normal": "JOB ENDED NORMALLY" in text,
        "hof": hof[-1] if hof else None,
        "gn": gn[-1] if gn else None,
        "ncyc": len(cyc),
        "g1": cyc[0][1] if cyc else None,
        "h1": cyc[0][2] if cyc else None,
        "bad": bad,
    }


def main(tags):
    print("%-22s | %-34s | %-34s" % ("case", "C++", "official"))
    print("-" * 100)
    for mol in MOLS:
        for kwtag, kw in KEYWORDS:
            tag = "%s_%s" % (mol, kwtag)
            if tags and tag not in tags:
                continue
            c = run(CPP, tag, kw, MOLS[mol])
            r = run(REF, tag, kw, MOLS[mol])

            def fmt(v):
                return ("exit=%s norm=%s HOF=%s GN=%s cyc=%s g1=%s"
                        % (v["exit"], int(v["normal"]), v["hof"], v["gn"],
                           v["ncyc"], v["g1"]))

            flag = ""
            if c["exit"] != 0 or not c["normal"]:
                flag = "  <<< CPP FAILED"
            elif c["hof"] != r["hof"]:
                flag = "  <<< HOF DIFFERS"
            print("%-22s | %-34s | %-34s%s" % (tag, fmt(c), fmt(r), flag))
            if c["bad"]:
                print("    cpp messages:", c["bad"])
            if r["bad"]:
                print("    ref messages:", r["bad"])
    return 0


if __name__ == "__main__":
    sys.exit(main([a for a in sys.argv[1:] if a != "--all"]))
