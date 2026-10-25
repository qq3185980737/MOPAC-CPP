#!/usr/bin/env python
"""run_cosmo_moz.py — MOZYME + COSMO regression against the official reference.

Runs a set of MOZYME+COSMO (EPS=78.4) decks with the C++ port and the official
Fortran reference and compares HOF / DIELECTRIC ENERGY / termination.

Usage: python run_cosmo_moz.py [--all]
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
ACETATE = """C -1.2000 0.0000 0.0000
C 0.1000 0.0000 0.0000
O 1.1000 0.6500 0.0000
O 1.1000 -0.6500 0.0000
H -1.6500 0.9500 0.0000
H -1.6500 -0.4800 0.8300
H -1.6500 -0.4800 -0.8300
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
CYCLOPENTADIENE = """C 0.0000 1.2080 0.0000
C 1.1489 0.3733 0.0000
C 0.7100 -0.9769 0.0000
C -0.7100 -0.9769 0.0000
C -1.1489 0.3733 0.0000
H 0.5500 2.0800 0.0000
H -0.5500 2.0800 0.0000
H 2.1800 0.7080 0.0000
H 1.3470 -1.8540 0.0000
H -1.3470 -1.8540 0.0000
H -2.1800 0.7080 0.0000
"""
NAPHTHALENE = """C 1.4000 0.0000 0.0000
C 0.7000 1.2124 0.0000
C -0.7000 1.2124 0.0000
C -1.4000 0.0000 0.0000
C -0.7000 -1.2124 0.0000
C 0.7000 -1.2124 0.0000
C 1.4000 2.4249 0.0000
C 0.7000 3.6373 0.0000
C -0.7000 3.6373 0.0000
C -1.4000 2.4249 0.0000
H 2.4840 0.0000 0.0000
H -2.4840 0.0000 0.0000
H -1.2420 -2.1510 0.0000
H 1.2420 -2.1510 0.0000
H 2.4840 2.4249 0.0000
H -2.4840 2.4249 0.0000
H 1.2420 4.5760 0.0000
H -1.2420 4.5760 0.0000
"""

# name: (keyword line, geometry)
CASES = {
    # MOZYME + COSMO 1SCF
    "water_moz_cosmo":    ("PM7 MOZYME 1SCF EPS=78.4", WATER),
    "meoh_moz_cosmo":     ("PM7 MOZYME 1SCF EPS=78.4", METHANOL),
    "acetate_moz_cosmo":  ("PM7 MOZYME 1SCF CHARGE=-1 EPS=78.4", ACETATE),
    "benzene_moz_cosmo":  ("PM7 MOZYME 1SCF EPS=78.4", BENZENE),
    "cpdiene_moz_cosmo":  ("PM7 MOZYME 1SCF EPS=78.4", CYCLOPENTADIENE),
    "naphth_moz_cosmo":   ("PM7 MOZYME 1SCF EPS=78.4", NAPHTHALENE),
    # ordinary (non-MOZYME) COSMO: 1SCF and optimisation
    "water_plain_cosmo":     ("PM7 1SCF EPS=78.4", WATER),
    "meoh_plain_cosmo":      ("PM7 1SCF EPS=78.4", METHANOL),
    "water_plain_cosmo_opt": ("PM7 EPS=78.4", WATER),
    "meoh_plain_cosmo_opt":  ("PM7 EPS=78.4", METHANOL),
    # MOZYME + COSMO optimisation (limited by the known MOZYME gradient issue)
    "water_moz_cosmo_opt": ("PM7 MOZYME EPS=78.4", WATER),
    "meoh_moz_cosmo_opt":  ("PM7 MOZYME EPS=78.4", METHANOL),
}

RE_HOF = re.compile(r"HEAT OF FORMATION\s*=\s*(-?\d+\.\d+)")
RE_ED = re.compile(r"DIELECTRIC ENERGY\s*=\s*(-?\d+\.\d+)")


def run(exe, tag, kw, geo):
    sub = os.path.join(ROOT, "cosmo_moz", tag, "cpp" if exe == CPP else "ref")
    os.makedirs(sub, exist_ok=True)
    with open(os.path.join(sub, tag + ".mop"), "w") as f:
        f.write(kw + "\n" + tag + "\n\n" + geo)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    with open(os.path.join(sub, "stdout.txt"), "wb") as fo, \
         open(os.path.join(sub, "stderr.txt"), "wb") as fe:
        p = subprocess.run([exe, tag + ".mop"], cwd=sub, env=env,
                           stdout=fo, stderr=fe, timeout=3600)
    text = ""
    out = os.path.join(sub, tag + ".out")
    if os.path.exists(out):
        with open(out, "r", errors="replace") as f:
            text = f.read()
    hof = RE_HOF.findall(text)
    ed = RE_ED.findall(text)
    return {
        "exit": p.returncode,
        "hof": hof[-1] if hof else None,
        "ediel": ed[-1] if ed else None,
        "normal": "JOB ENDED NORMALLY" in text,
    }


def main(tags):
    ok_all = True
    for tag in tags:
        kw, geo = CASES[tag]
        c = run(CPP, tag, kw, geo)
        r = run(REF, tag, kw, geo)
        ok = (c["hof"] == r["hof"] and c["ediel"] == r["ediel"]
              and c["exit"] == r["exit"] == 0 and c["normal"] and r["normal"])
        ok_all = ok_all and ok
        print("%-22s cpp HOF=%-12s EDIEL=%-10s | ref HOF=%-12s EDIEL=%-10s | %s"
              % (tag, c["hof"], c["ediel"], r["hof"], r["ediel"],
                 "PASS" if ok else "FAIL"))
        if not ok:
            print("   cpp:", c)
            print("   ref:", r)
    print("RESULT:", "ALL PASS" if ok_all else "FAILURES")
    return 0 if ok_all else 1


if __name__ == "__main__":
    tags = list(CASES)
    if len(sys.argv) > 1 and sys.argv[1] != "--all":
        tags = [a for a in sys.argv[1:] if a in CASES]
    sys.exit(main(tags))
