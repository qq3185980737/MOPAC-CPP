// moz_cphf.h — first building blocks of the exact analytic MOZYME gradient
// (CPHF/LMO-response) project.  Diagnostics only; opt-in via MOPAC_MOZ_WDBG.
#pragma once
#include <vector>

// Validates the dense-LMO reconstruction against the packed MOZYME density and
// reports Tr(W), max|W-W^T| for W = C F C^T (F = packed AO Fock).
void moz_cphf_selftest(const std::vector<double>& packed_fock);
