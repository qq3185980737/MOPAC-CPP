# Performance benchmark decks (C40H82, C80H162, C120H242)

Single-point PM7 RHF workloads used for the performance results of the SPE
manuscript (Section 3.6). All geometries are physically valid all-trans
alkanes (C-C-C = 109.47 degrees, C-C = 1.532 A, C-H = 1.093 A).

## Files

- `c40h82_pm7_1scf.mop`  - C40H82, 122 atoms
- `c80h162_pm7_1scf.mop` - C80H162, 242 atoms
- `c120h242_pm7_1scf.mop` - C120H242, 362 atoms (722 basis functions)

## Protocol

Single thread, MKL single-threaded, warm-up followed by five timed runs
(median reported on the paper):

```
set OMP_NUM_THREADS=1
set MKL_NUM_THREADS=1
mopac.exe c80h162_pm7_1scf.mop
```

For the deterministic-parallel (FBDR) runs set `MOPAC_DSPAR=1`
(optionally `MOPAC_DSPAR_NB=8`) and `OMP_NUM_THREADS=8`; the current build
needs `MOPAC_DSPAR_SCF=1` for the SCF path as well. The FINAL HEAT must be
bitwise identical for every thread count.

## Reference results (2026-11, AMD Ryzen 7 7840HS, Windows x64)

Port versus the official MOPAC2016 rebuild (ifort + MKL), single thread:

| Workload | Official rebuild | C++ port | Ratio |
|---|---|---|---|
| C80H162 | 4.30 s | 2.05 s | 2.1x |
| C120H242 | 9.34 s | 4.29 s | 2.18x |

FBDR thread scaling on the current AVX-512/OpenMP build is
performance-neutral (0.8-1.0x at 4-16 threads); its value is determinism,
not speed.
