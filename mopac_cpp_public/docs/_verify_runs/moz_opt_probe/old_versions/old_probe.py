#!/usr/bin/env python
"""old_probe.py — run the MOZYME-optimisation probe with an alternative C++
executable (e.g. the old GitHub version) and compare with the official reference.

Usage:
    python old_probe.py <cpp_exe> [tag_subset...]
"""
import os
import re
import subprocess
import sys
import importlib.util

ROOT = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(ROOT, "bin", "mopac_ref.exe")

spec = importlib.util.spec_from_file_location(
    "probe", os.path.join(ROOT, "run_moz_opt_probe.py"))
probe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(probe)

RE_HOF = probe.RE_HOF
RE_GN = probe.RE_GN
RE_CYC = probe.RE_CYC


def run(exe, tag, kw, geo, outroot):
    sub = os.path.join(outroot, tag, "cpp" if exe != REF else "ref")
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
    return {"exit": rc, "normal": "JOB ENDED NORMALLY" in text,
            "hof": hof[-1] if hof else None,
            "gn": gn[-1] if gn else None,
            "ncyc": len(cyc), "g1": cyc[0][1] if cyc else None}


def main():
    cpp = sys.argv[1]
    tags = sys.argv[2:] or ["water_default", "water_COSMO",
                            "benzene_default", "benzene_COSMO"]
    outroot = os.path.join(ROOT, "moz_opt_old")
    print("cpp exe:", cpp)
    print("%-20s | %-30s | %-30s" % ("case", "old C++ / current C++", "official"))
    print("-" * 90)
    for tag in tags:
        mol, kwtag = tag.split("_", 1)
        if mol not in probe.MOLS:
            continue
        kw = dict(probe.KEYWORDS)[kwtag]
        c = run(cpp, tag, kw, probe.MOLS[mol], outroot)
        r = run(REF, tag, kw, probe.MOLS[mol], outroot)
        print("%-20s | g1=%-9s HOF=%-11s GN=%-8s cyc=%s | g1=%-8s HOF=%-11s GN=%s cyc=%s"
              % (tag, c["g1"], c["hof"], c["gn"], c["ncyc"],
                 r["g1"], r["hof"], r["gn"], r["ncyc"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
