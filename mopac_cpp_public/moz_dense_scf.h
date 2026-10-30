// moz_dense_scf.h — internal reimplementation of the MOZYME SCF iteration.
//
// Instead of the incremental/DIIS/level-shift machinery of iter_for_MOZYME, this
// module runs a variational LMO-SCF on a dense occupied coefficient matrix C:
//   P = 2 C C^T  ->  packed density  ->  Fock  ->  residual R = F C - C (C^T F C)
//   C <- orthonormalised(C - eta R), with backtracking on the energy
// until max|R| < tol.  C is kept in static storage across calls, so the solution
// is continued from the previous geometry (branch memory).
//
// Opt-in: MOPAC_MOZ_DENSE_SCF=1 (MOPAC_MOZ_DENSE_MAXIT, MOPAC_MOZ_DENSE_TOL).
#pragma once

// Returns true if the dense SCF handled the call; ee_out = electronic energy
// (the quantity iter_for_MOZYME normally returns through its first argument).
bool moz_dense_scf(double& ee_out);
