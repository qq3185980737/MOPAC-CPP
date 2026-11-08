# MOPAC-CPP v1.2.0 release notes

Date: 2026-11-09

## Scope

This release freezes the state of the C++17 translation used for the SPE
manuscript, including the opt-in UHF/MD fixes and the updated verification
documentation. The default execution path remains the bitwise-comparable one.

## Changes since v1.1.0 (highlights)

- UHF: optional defect fixes behind `MOPAC_UHF_FIX` (beta residual in the SCF
  exit test, beta level-shift gap from the beta spectrum, FULLDIAG-aware beta
  diagonalization, open-shell SCFCRT tightening). Default OFF.
- UHF: `MOPAC_DUMP_UHF` diagnostic dump for offline per-spin analysis.
- MD: high-precision energy print (`MOPAC_MD_HP`), dcart gradient dump
  (`MOPAC_GRADDUMP`), and delta-ML hook (`MOPAC_DML_MODEL`, experimental).
- MOZYME: dense-SCF prototypes and negative-result probes (not on the default
  path); see `experiments/` locally and the analysis notes.
- Benchmarks: C40/C80/C120 PM7 1SCF performance decks and timing protocol.
- Docs: updated `FORMULAS_AND_VERIFICATION.docx`, added
  `MOZYME_MD_conservation.md`, `.gitignore`, `CITATION.cff`.

## Verification status

- Build: MSVC 2022 x64 Release, Intel MKL 2025.3, SSE2 baseline.
- Released regression corpus: 327 paired inputs, **320 bitwise identical**,
  7 EF-optimization 1-ULP endpoint mismatches (documented build sensitivity),
  22 MD/DRC inputs without a Fortran `.arc` counterpart.
- Smoke test: `c_hatom_plain.mop` runs to normal termination with
  `FINAL HEAT OF FORMATION = 52.10200 kcal/mol`.

## Known limitations

- Cross-platform bitwise equality is not claimed (single Windows x64 build pair).
- The UHF/MD/delta-ML paths are experimental and off by default.
- Zenodo DOI is assigned after deposit; update the manuscript's data-availability
  statement with the DOI once available.
