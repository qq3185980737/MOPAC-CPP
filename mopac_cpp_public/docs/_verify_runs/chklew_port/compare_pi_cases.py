#!/usr/bin/env python
"""compare_pi_cases.py — fidelity check of the rewritten CHKLEW pi framework.

Runs a set of MOZYME single-point decks with the C++ port and the official
Fortran reference and reports whether the heat of formation agrees exactly.
This exercises the ene/yne (open-ended), aromatic six-ring, fused-ring and
five-membered-ring branches of chklew.F90 289-618.

S-containing systems are deliberately absent: MOZYME + sulphur is broken in the
current C++ port for reasons unrelated to CHKLEW (thiophene/methanethiol give
wild energies with both the old and the new CHKLEW, while the plain PM7 path
matches the reference exactly; CHKLEW's GEOCHK dumps for thiophene are already
identical to the reference).  See docs/TRANSLATION_GAPS.md §95 for details.

Usage: python compare_pi_cases.py
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(ROOT, "bin")
CPP = os.path.join(BIN, "mopac_cpp.exe")
REF = os.path.join(BIN, "mopac_ref.exe")
CASES = os.path.join(ROOT, "pi_cases")

DECKS = {
    # open-ended (ene) pi systems
    "ethene": """PM7 MOZYME 1SCF
ethene mozyme

C 0.0000 0.0000 0.6670
C 0.0000 0.0000 -0.6670
H 0.0000 0.9230 1.2380
H 0.0000 -0.9230 1.2380
H 0.0000 0.9230 -1.2380
H 0.0000 -0.9230 -1.2380
""",
    # yne: acetylenic double pi bond (F90 495-497)
    "acetylene": """PM7 MOZYME 1SCF
acetylene mozyme

C 0.0000 0.0000 0.6000
C 0.0000 0.0000 -0.6000
H 0.0000 0.0000 1.6600
H 0.0000 0.0000 -1.6600
""",
    # conjugated ene chain
    "butadiene": """PM7 MOZYME 1SCF
trans-butadiene mozyme

C 0.0000 0.0000 0.0000
C 1.3400 0.0000 0.0000
C 2.0900 1.2900 0.0000
C 3.4300 1.2900 0.0000
H -0.5000 0.9300 0.0000
H -0.5000 -0.9300 0.0000
H 1.8800 -0.9300 0.0000
H 1.5700 2.2400 0.0000
H 2.8900 2.2300 0.0000
H 2.8900 0.3600 0.0000
""",
    # heteroatom aromatic ring
    "pyridine": """PM7 MOZYME 1SCF
pyridine mozyme

N 1.3970 0.0000 0.0000
C 0.6985 1.2100 0.0000
C -0.6985 1.2100 0.0000
C -1.3970 0.0000 0.0000
C -0.6985 -1.2100 0.0000
C 0.6985 -1.2100 0.0000
H 1.2405 2.1490 0.0000
H -1.2405 2.1490 0.0000
H -2.4810 0.0000 0.0000
H -1.2405 -2.1490 0.0000
H 1.2405 -2.1490 0.0000
""",
    # fused six-membered rings (ring phase, 1040 restart)
    "naphthalene": """PM7 MOZYME 1SCF
naphthalene mozyme

C 1.4000 0.0000 0.0000
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
""",
    # five-membered ring (ring5 + F90 252-288)
    "cyclopentadiene": """PM7 MOZYME 1SCF
cyclopentadiene mozyme

C 0.0000 1.2080 0.0000
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
""",
    # fused five- plus six-membered ring (F90 252-288 opening)
    "indole": """PM7 MOZYME 1SCF
indole mozyme

C 1.4000 0.0000 0.0000
C 0.7000 1.2124 0.0000
C -0.7000 1.2124 0.0000
C -1.4000 0.0000 0.0000
C -0.7000 -1.2124 0.0000
C 0.7000 -1.2124 0.0000
C 1.1320 2.5434 0.0000
N 0.0000 3.3664 0.0000
C -1.1320 2.5434 0.0000
H 2.4840 0.0000 0.0000
H -2.4840 0.0000 0.0000
H -1.2420 -2.1510 0.0000
H 1.2420 -2.1510 0.0000
H 1.9800 2.8500 0.0000
H -1.9800 2.8500 0.0000
H 0.0000 4.3800 0.0000
""",
    # three fused six-membered rings (ring phase, repeated 1040 restart)
    "anthracene": """PM7 MOZYME 1SCF
anthracene mozyme

C -1.2125 0.7000 0.0000
C -2.4249 1.4000 0.0000
C -3.6373 0.7000 0.0000
C -3.6373 -0.7000 0.0000
C -2.4249 -1.4000 0.0000
C -1.2125 -0.7000 0.0000
C 0.0000 1.4000 0.0000
C 1.2125 0.7000 0.0000
C 1.2125 -0.7000 0.0000
C 0.0000 -1.4000 0.0000
C 2.4249 1.4000 0.0000
C 3.6373 0.7000 0.0000
C 3.6373 -0.7000 0.0000
C 2.4249 -1.4000 0.0000
H -2.4249 2.4800 0.0000
H -4.5760 1.2420 0.0000
H -4.5760 -1.2420 0.0000
H -2.4249 -2.4800 0.0000
H 0.0000 2.4800 0.0000
H 0.0000 -2.4800 0.0000
H 2.4249 2.4800 0.0000
H 4.5760 1.2420 0.0000
H 4.5760 -1.2420 0.0000
H 2.4249 -2.4800 0.0000
""",
}


def write_decks():
    os.makedirs(CASES, exist_ok=True)
    for tag, text in DECKS.items():
        with open(os.path.join(CASES, tag + ".mop"), "w", newline="\n") as f:
            f.write(text)


def run(exe, tag):
    sub = os.path.join(CASES, tag, "cpp" if exe == CPP else "ref")
    os.makedirs(sub, exist_ok=True)
    with open(os.path.join(sub, tag + ".mop"), "w") as f:
        f.write(DECKS[tag])
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    with open(os.path.join(sub, "stdout.txt"), "wb") as fo, \
         open(os.path.join(sub, "stderr.txt"), "wb") as fe:
        p = subprocess.run([exe, tag + ".mop"], cwd=sub, env=env,
                           stdout=fo, stderr=fe, timeout=1800)
    out = os.path.join(sub, tag + ".out")
    text = ""
    if os.path.exists(out):
        with open(out, "r", errors="replace") as f:
            text = f.read()
    hof = None
    import re
    m = re.findall(r"HEAT OF FORMATION\s*=\s*(-?\d+\.\d+)", text)
    if m:
        hof = m[-1]
    ok = "JOB ENDED NORMALLY" in text
    return p.returncode, hof, ok


def main():
    write_decks()
    all_ok = True
    for tag in DECKS:
        ec, hc, okc = run(CPP, tag)
        er, hr, okr = run(REF, tag)
        same = (hc == hr)
        all_ok = all_ok and same and ec == er == 0 and okc and okr
        print("{0:16s} cpp={1:>12} ref={2:>12} exit={3}/{4} normal={5}/{6} {7}"
              .format(tag, hc, hr, ec, er, okc, okr, "OK" if same else "DIFF"))
    print("RESULT:", "ALL IDENTICAL" if all_ok else "DIFFERENCES PRESENT")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
