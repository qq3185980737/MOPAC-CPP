#!/usr/bin/env python
"""run_gate.py — run MOPAC decks with the C++ port and the Fortran reference,
then compare the gate metrics against the recorded baseline.

Usage:
    python run_gate.py deck1.mop [deck2.mop ...]
    python run_gate.py --all

For each deck (tag = stem), the two binaries are run in
    <root>/run/<tag>/cpp/<tag>.mop   -> cpp.out
    <root>/run/<tag>/ref/<tag>.mop   -> ref.out
The output files are kept for inspection, and one verdict line per deck is
printed.  A deck PASSES when the C++ HOF (and ZPE, if the baseline records one)
equals the value in baseline.json.  Exact agreement with the Fortran reference
is reported too, but the reference ZPE is allowed to differ where the baseline
records a known dcart/THEMO precision divergence.
"""
import json
import os
import re
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(ROOT, "bin")
CPP = os.path.join(BIN, "mopac_cpp.exe")
REF = os.path.join(BIN, "mopac_ref.exe")

RE_HOF = re.compile(r"HEAT OF FORMATION\s*=\s*(-?\d+\.\d+)")
RE_ZPE = re.compile(r"ZERO POINT ENERGY\s+(-?\d+\.\d+)")

ALL_DECKS = [
    "water_bonds", "benzene_pi", "water_thermo",
    "water_force", "water_mozyme", "benzene_mozyme",
    # chklew pi-framework fidelity cases (pi_cases/):
    "ethene", "acetylene", "butadiene", "pyridine", "naphthalene",
    "cyclopentadiene", "indole", "anthracene",
]


def run_one(exe, tag, mop, sub):
    os.makedirs(sub, exist_ok=True)
    shutil.copyfile(mop, os.path.join(sub, tag + ".mop"))
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    out_path = os.path.join(sub, tag + ".out")
    with open(os.path.join(sub, "stdout.txt"), "wb") as fout, \
         open(os.path.join(sub, "stderr.txt"), "wb") as ferr:
        proc = subprocess.run([exe, tag + ".mop"], cwd=sub, env=env,
                              stdout=fout, stderr=ferr, timeout=3600)
    text = ""
    if os.path.exists(out_path):
        with open(out_path, "r", errors="replace") as f:
            text = f.read()
    hof = RE_HOF.findall(text)
    zpe = RE_ZPE.findall(text)
    return {
        "exit": proc.returncode,
        "hof": hof[-1] if hof else None,
        "zpe": zpe[-1] if zpe else None,
        "out": out_path,
    }


def main(tags):
    with open(os.path.join(ROOT, "baseline.json")) as f:
        baseline = json.load(f)
    ok_all = True
    results = {}
    for tag in tags:
        mop = None
        for d in ("decks", "pi_cases", "run"):
            cand = os.path.join(ROOT, d, tag + ".mop")
            if os.path.exists(cand):
                mop = cand
                break
        if mop is None:
            print("[{0}] SKIP (missing deck {1})".format(tag, tag + ".mop"))
            ok_all = False
            continue
        base = os.path.join(ROOT, "run", tag)
        c = run_one(CPP, tag, mop, os.path.join(base, "cpp"))
        r = run_one(REF, tag, mop, os.path.join(base, "ref"))
        exp = baseline.get(tag, {})
        ok_hof = c["hof"] == exp.get("hof")
        ok_zpe = c["zpe"] == exp.get("zpe")
        ref_hof_match = c["hof"] == r["hof"]
        ref_zpe_match = c["zpe"] == r["zpe"]
        ok = ok_hof and ok_zpe and c["exit"] == 0
        ok_all = ok_all and ok
        print("[{0}] cpp HOF={1} ZPE={2} | ref HOF={3} ZPE={4} | "
              "baseline={5} | {6}".format(
                  tag, c["hof"], c["zpe"], r["hof"], r["zpe"],
                  "OK" if ok else "CHANGED", "PASS" if ok else "FAIL"))
        if not ok_hof:
            print("    HOF changed: baseline={0!r} now={1!r} (ref={2!r})"
                  .format(exp.get("hof"), c["hof"], r["hof"]))
        if not ok_zpe:
            print("    ZPE changed: baseline={0!r} now={1!r} (ref={2!r})"
                  .format(exp.get("zpe"), c["zpe"], r["zpe"]))
        if c["exit"] != 0:
            print("    cpp exit code = {0}".format(c["exit"]))
        results[tag] = {"cpp": c, "ref": r, "baseline_ok": ok,
                        "ref_hof_match": ref_hof_match,
                        "ref_zpe_match": ref_zpe_match}
    with open(os.path.join(ROOT, "run", "last_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("RESULT:", "ALL PASS" if ok_all else "FAILURES PRESENT")
    return 0 if ok_all else 1


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        sys.exit(2)
    if args == ["--all"]:
        tags = ALL_DECKS
    else:
        tags = [os.path.splitext(os.path.basename(a))[0] for a in args]
    sys.exit(main(tags))
