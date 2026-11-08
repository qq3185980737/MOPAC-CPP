# Verification summary — MOPAC-CPP v1.2.0 (2026-11-09)

- Commit: `8dc5149` (main), tag `v1.2.0`
- Build: MSVC 2022 x64 Release, Intel MKL 2025.3, SSE2 baseline.
  - `MSBuild MOPAC_Cpp.sln /p:Configuration=Release /p:Platform=x64` -> **succeeded**
  - Output: `mopac_cpp_public/x64/Release/mopac.exe`
- Smoke test: `c_hatom_plain.mop` (single point)
  - `FINAL HEAT OF FORMATION = 52.10200 KCAL/MOL = 217.99477 KJ/MOL`
  - `MOPAC Job ended normally`
- Released regression corpus (`mopac_cpp_public/regression`, `cmp_corpus.py`):
  - 327 paired inputs, **320 bitwise identical**
  - 7 EF-optimization endpoint mismatches at 1-ULP level (documented build sensitivity)
  - 22 MD/DRC inputs without Fortran `.arc` counterpart
  - `missing_cxx_output: 0`
- Scope: default single-threaded bitwise path unchanged; UHF/MD/delta-ML paths are
  opt-in and off by default.
- Next: GitHub Release for `v1.2.0`, Zenodo deposit, then replace the DOI placeholder
  in the SPE manuscript data-availability statement.
