// moz_cphf.h — exact analytic MOZYME gradient project (CPHF / LMO response).
#pragma once
#include <vector>

// Validates the dense-LMO reconstruction against the packed MOZYME density and
// reports the Fock-matrix stationarity residual (max|F_ov|) and idempotency.
void moz_cphf_selftest(const std::vector<double>& packed_fock);

// Matrix-free CPHF/Z-vector core.  Requires a converged MOZYME SCF state.
// Builds F_ov = C_occ^T F C_vir, defines the response operator
//     A x = (eps_a - eps_i) x + [F(P + eps*dP(x)) - F(P)]_ov / eps
// (dP(x) = C_occ x C_vir^T + h.c., the Fock response being obtained with
// buildf on the perturbed density) and solves A z = -2 F_ov by conjugate
// gradients.  Opt-in via MOPAC_MOZ_CPHF=1; prints CPHF diagnostics.
void moz_cphf_run();
