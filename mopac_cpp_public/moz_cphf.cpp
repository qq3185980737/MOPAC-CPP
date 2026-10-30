// ============================================================================
// moz_cphf.cpp — exact analytic MOZYME gradient project (CPHF/LMO response).
//
// Stage 1: reconstruct dense occupied LMO matrix C and validate P = 2 C C^T
//          against the packed MOZYME density.
// Stage 2: S = I in NDDO, so there is no dS/dR (Pulay-overlap) term.
// Stage 3: the remaining gradient error is the density-response term
//          Tr((dE/dP) . dP/dR); here we quantify its starting point: the SCF
//          stationarity residual, i.e. the occupied-virtual block of the Fock
//          matrix in the LMO basis.
//
// Diagnostics only; opt-in via MOPAC_MOZ_WDBG.
// ============================================================================

#include "moz_cphf.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "ijbo.h"
#include "buildf.h"
#include <chrono>

namespace MOZYME_C {
extern int noccupied;
extern int nvirtual;
}

namespace {

// Build a dense LMO coefficient block (norbs x nmo) from the MOZYME storage.
// The occupied and virtual sets use parallel bookkeeping arrays:
//   occupied: icocc[nncf[i]+1 .. +ncf[i]], coefficients cocc[ncocc[i] + ...]
//   virtual : icvir[nnce[i]+1 .. +nce[i]], coefficients cvir[ncvir[i] + ...]
void build_dense_lmo(const std::vector<int>& iatom,
                     const std::vector<int>& natom,
                     const std::vector<int>& start,
                     const std::vector<double>& coef,
                     const std::vector<int>& first_list,
                     const std::vector<int>& count_list,
                     int nmo, std::vector<double>& C) {
  using namespace MOZYME_C;
  const int nr = molkst_C::norbs;
  C.assign((size_t)nr * nmo, 0.0);
  for (int i = 1; i <= nmo; ++i) {
    const int base = start[i];
    int ja = 0;
    for (int jj = 1; jj <= count_list[i]; ++jj) {
      const int atom = iatom[first_list[i] + jj];
      const int nj = iorbs[atom];
      const int first_ao = common_arrays_C::nfirst[atom];
      for (int mu = 1; mu <= nj; ++mu) {
        const int row = first_ao - 1 + mu - 1;
        C[(size_t)row * nmo + (i - 1)] = coef[base + ja + mu];
      }
      ja += nj;
    }
    (void)natom;
  }
}

int ijbo_safe(int i, int j) { return ijbo(i, j); }

// packed (MOZYME) -> dense AO matrix; see helecz/dcart conventions
void dense_from_packed(const std::vector<double>& packed,
                       std::vector<double>& dense) {
  using namespace MOZYME_C;
  const int nr = molkst_C::norbs;
  dense.assign((size_t)nr * nr, 0.0);
  for (int a = 1; a <= molkst_C::numat; ++a) {
    const int fa = common_arrays_C::nfirst[a];
    for (int b = 1; b <= a; ++b) {
      const int k = ijbo_safe(a, b);
      if (k < 0) continue;
      const int fb = common_arrays_C::nfirst[b];
      if (a == b) {
        int idx = k;
        for (int mu = 1; mu <= iorbs[a]; ++mu)
          for (int nu = 1; nu <= mu; ++nu) {
            ++idx;
            const double v = packed[idx];
            dense[(size_t)(fa - 1 + mu - 1) * nr + (fb - 1 + nu - 1)] = v;
            dense[(size_t)(fb - 1 + nu - 1) * nr + (fa - 1 + mu - 1)] = v;
          }
      } else {
        int idx = k;
        for (int mu = 1; mu <= iorbs[a]; ++mu)
          for (int nu = 1; nu <= iorbs[b]; ++nu) {
            ++idx;
            const double v = packed[idx];
            dense[(size_t)(fa - 1 + mu - 1) * nr + (fb - 1 + nu - 1)] = v;
            dense[(size_t)(fb - 1 + nu - 1) * nr + (fa - 1 + mu - 1)] = v;
          }
      }
    }
  }
}

}  // namespace

void moz_cphf_selftest(const std::vector<double>& packed_fock) {
  using namespace MOZYME_C;
  const int nr = molkst_C::norbs;
  const int nocc = MOZYME_C::noccupied;
  const int nvir = MOZYME_C::nvirtual;
  if (nr <= 0 || nocc <= 0) {
    std::fprintf(stderr, "MOZCPHF norbs/nocc invalid (%d,%d)\n", nr, nocc);
    return;
  }

  std::vector<double> Cocc;
  build_dense_lmo(icocc, ncf, ncocc, cocc, nncf, ncf, nocc, Cocc);

  // ---- stage 1: P = 2 C C^T vs packed density ----
  std::vector<double> P((size_t)nr * nr, 0.0);
  for (int mu = 0; mu < nr; ++mu)
    for (int nu = 0; nu < nr; ++nu) {
      double s = 0.0;
      for (int i = 0; i < nocc; ++i) s += Cocc[(size_t)mu * nocc + i] * Cocc[(size_t)nu * nocc + i];
      P[(size_t)mu * nr + nu] = 2.0 * s;
    }
  std::vector<double> p_dense;
  dense_from_packed(common_arrays_C::p, p_dense);
  double maxerr = 0.0;
  for (size_t i = 0; i < P.size(); ++i)
    maxerr = std::max(maxerr, std::fabs(P[i] - p_dense[i]));

  // ---- idempotency of the single-particle density ----
  double idem = 0.0;
  for (int mu = 0; mu < nr; ++mu)
    for (int nu = 0; nu < nr; ++nu) {
      double s = 0.0;
      for (int k = 0; k < nr; ++k) s += P[(size_t)mu * nr + k] * P[(size_t)k * nr + nu];
      idem = std::max(idem, std::fabs(s - 2.0 * P[(size_t)mu * nr + nu]));
    }

  std::fprintf(stderr, "MOZCPHF norb=%d nocc=%d nvir=%d\n", nr, nocc, nvir);
  std::fprintf(stderr, "MOZCPHF density max|P_dense-P_packed| = %.3e  max|P^2-2P| = %.3e\n",
               maxerr, idem);

  // ---- stage 3: Fock in the LMO basis; occ-virt stationarity residual ----
  std::vector<double> F;
  dense_from_packed(packed_fock, F);
  auto mo_block = [&](const std::vector<double>& C1, int n1,
                      const std::vector<double>& C2, int n2) {
    // returns C1^T F C2 (n1 x n2)
    std::vector<double> out((size_t)n1 * n2, 0.0);
    std::vector<double> FC2((size_t)nr * n2, 0.0);
    for (int mu = 0; mu < nr; ++mu)
      for (int j = 0; j < n2; ++j) {
        double s = 0.0;
        for (int nu = 0; nu < nr; ++nu) s += F[(size_t)mu * nr + nu] * C2[(size_t)nu * n2 + j];
        FC2[(size_t)mu * n2 + j] = s;
      }
    for (int i = 0; i < n1; ++i)
      for (int j = 0; j < n2; ++j) {
        double s = 0.0;
        for (int mu = 0; mu < nr; ++mu) s += C1[(size_t)mu * n1 + i] * FC2[(size_t)mu * n2 + j];
        out[(size_t)i * n2 + j] = s;
      }
    return out;
  };

  std::vector<double> Foo = mo_block(Cocc, nocc, Cocc, nocc);
  double foo_max = 0.0, foo_diag = 0.0;
  for (int i = 0; i < nocc; ++i) {
    foo_diag += Foo[(size_t)i * nocc + i];
    for (int j = 0; j < nocc; ++j) foo_max = std::max(foo_max, std::fabs(Foo[(size_t)i * nocc + j]));
  }
  std::fprintf(stderr, "MOZCPHF Tr(F_oo)=%.6f  max|F_oo|=%.6f\n", foo_diag, foo_max);

  if (nvir > 0) {
    std::vector<double> Cvir;
    build_dense_lmo(icvir, nce, ncvir, cvir, nnce, nce, nvir, Cvir);
    std::vector<double> Fov = mo_block(Cocc, nocc, Cvir, nvir);
    std::vector<double> Fvv = mo_block(Cvir, nvir, Cvir, nvir);
    double fov_max = 0.0, fov_rms = 0.0;
    for (double v : Fov) {
      fov_max = std::max(fov_max, std::fabs(v));
      fov_rms += v * v;
    }
    fov_rms = std::sqrt(fov_rms / (double)Fov.size());
    double fvv_max = 0.0;
    for (double v : Fvv) fvv_max = std::max(fvv_max, std::fabs(v));
    std::fprintf(stderr,
                 "MOZCPHF stationarity: max|F_ov|=%.6f rms|F_ov|=%.6f "
                 "(max|F_oo|=%.6f max|F_vv|=%.6f)\n",
                 fov_max, fov_rms, foo_max, fvv_max);
  }
}

// ============================================================================
// Matrix-free CPHF/Z-vector core (appended): moz_cphf_run().
// Requires a converged MOZYME SCF state (globals p, partf, cocc/cvir storage).
// ============================================================================

#include "buildf.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

// dense (norbs x norbs) -> packed MOZYME layout (same blocks as p / f)
void packed_from_dense(const std::vector<double>& dense,
                       std::vector<double>& packed) {
  using namespace MOZYME_C;
  const int nr = molkst_C::norbs;
  packed.assign(common_arrays_C::p.size(), 0.0);
  for (int a = 1; a <= molkst_C::numat; ++a) {
    const int fa = common_arrays_C::nfirst[a];
    for (int b = 1; b <= a; ++b) {
      const int k = ijbo(a, b);
      if (k < 0) continue;
      const int fb = common_arrays_C::nfirst[b];
      if (a == b) {
        int idx = k;
        for (int mu = 1; mu <= iorbs[a]; ++mu)
          for (int nu = 1; nu <= mu; ++nu) {
            ++idx;
            packed[idx] = dense[(size_t)(fa - 1 + mu - 1) * nr + (fb - 1 + nu - 1)];
          }
      } else {
        int idx = k;
        for (int mu = 1; mu <= iorbs[a]; ++mu)
          for (int nu = 1; nu <= iorbs[b]; ++nu) {
            ++idx;
            packed[idx] = dense[(size_t)(fa - 1 + mu - 1) * nr + (fb - 1 + nu - 1)];
          }
      }
    }
  }
}

}  // namespace

void moz_cphf_run() {
  using namespace MOZYME_C;
  const int nr = molkst_C::norbs;
  const int nocc = MOZYME_C::noccupied;
  const int nvir = MOZYME_C::nvirtual;
  if (nr <= 0 || nocc <= 0 || nvir <= 0) return;
  const size_t n = (size_t)nocc * nvir;

  // ---- dense LMO blocks ----
  std::vector<double> Cocc, Cvir;
  // reuse the stage-1 builders (same TU, anonymous namespace)
  build_dense_lmo(icocc, ncf, ncocc, cocc, nncf, ncf, nocc, Cocc);
  build_dense_lmo(icvir, nce, ncvir, cvir, nnce, nce, nvir, Cvir);

  // ---- converged Fock and its occ-virt block ----
  std::vector<double> fcur(common_arrays_C::p.size(), 0.0);
  buildf(fcur, MOZYME_C::partf, 0);
  std::vector<double> Fd;
  dense_from_packed(fcur, Fd);

  auto project = [&](const std::vector<double>& F, const std::vector<double>& C1,
                     int n1, const std::vector<double>& C2, int n2,
                     std::vector<double>& M) {
    std::vector<double> FC2((size_t)nr * n2, 0.0);
    for (int mu = 0; mu < nr; ++mu)
      for (int j = 0; j < n2; ++j) {
        double s = 0.0;
        for (int nu = 0; nu < nr; ++nu) s += F[(size_t)mu * nr + nu] * C2[(size_t)nu * n2 + j];
        FC2[(size_t)mu * n2 + j] = s;
      }
    M.assign((size_t)n1 * n2, 0.0);
    for (int i = 0; i < n1; ++i)
      for (int j = 0; j < n2; ++j) {
        double s = 0.0;
        for (int mu = 0; mu < nr; ++mu) s += C1[(size_t)mu * n1 + i] * FC2[(size_t)mu * n2 + j];
        M[(size_t)i * n2 + j] = s;
      }
  };

  std::vector<double> Fov, Foo, Fvv;
  project(Fd, Cocc, nocc, Cvir, nvir, Fov);
  project(Fd, Cocc, nocc, Cocc, nocc, Foo);
  project(Fd, Cvir, nvir, Cvir, nvir, Fvv);

  double fov_rms = 0.0;
  for (double v : Fov) fov_rms += v * v;
  fov_rms = std::sqrt(fov_rms / (double)n);

  // diagonal orbital-energy model from the MO-basis Fock blocks
  std::vector<double> eo(nocc), ev(nvir);
  for (int i = 0; i < nocc; ++i) eo[i] = Foo[(size_t)i * nocc + i];
  for (int i = 0; i < nvir; ++i) ev[i] = Fvv[(size_t)i * nvir + i];

  const auto t0 = std::chrono::steady_clock::now();
  const std::vector<double> psave = common_arrays_C::p;
  const double epsP = 1.0e-2;

  auto apply_A = [&](const std::vector<double>& x, std::vector<double>& Ax) {
    // dP_dense = Cocc x Cvir^T (+ transpose)
    std::vector<double> dP((size_t)nr * nr, 0.0);
    for (int i = 0; i < nocc; ++i)
      for (int j = 0; j < nvir; ++j) {
        const double v = x[(size_t)i * nvir + j];
        if (v == 0.0) continue;
        for (int mu = 0; mu < nr; ++mu) {
          const double cm = Cocc[(size_t)mu * nocc + i] * v;
          if (cm == 0.0) continue;
          for (int nu = 0; nu < nr; ++nu) {
            const double cn = Cvir[(size_t)nu * nvir + j];
            if (cn != 0.0) dP[(size_t)mu * nr + nu] += cm * cn;
          }
        }
      }
    std::vector<double> dPpacked;
    packed_from_dense(dP, dPpacked);
    common_arrays_C::p = psave;
    for (size_t idx = 1; idx < dPpacked.size(); ++idx)
      common_arrays_C::p[idx] += epsP * dPpacked[idx];
    std::vector<double> f2(psave.size(), 0.0);
    buildf(f2, MOZYME_C::partf, 0);
    common_arrays_C::p = psave;
    std::vector<double> F2d;
    dense_from_packed(f2, F2d);
    std::vector<double> F2ov;
    project(F2d, Cocc, nocc, Cvir, nvir, F2ov);
    Ax.resize(n);
    for (int i = 0; i < nocc; ++i)
      for (int j = 0; j < nvir; ++j) {
        const size_t k = (size_t)i * nvir + j;
        Ax[k] = (ev[j] - eo[i]) * x[k] + (F2ov[k] - Fov[k]) / epsP;
      }
  };

  auto dot = [](const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
  };

  // RHS b = -2 F_ov
  std::vector<double> b(n), x(n, 0.0), r(n), pvec(n), Ax(n);
  for (size_t k = 0; k < n; ++k) b[k] = -2.0 * Fov[k];
  r = b;
  pvec = r;
  const double bnorm = std::sqrt(dot(b, b));
  double rsold = dot(r, r);
  int it = 0;
  const int maxit = (int)std::getenv("MOPAC_MOZ_CPHF_MAXIT") ?
      std::atoi(std::getenv("MOPAC_MOZ_CPHF_MAXIT")) : 100;
  const double tol = bnorm > 0.0 ? 1.0e-6 * bnorm : 0.0;
  for (it = 1; it <= maxit; ++it) {
    apply_A(pvec, Ax);
    const double denom = dot(pvec, Ax);
    if (std::fabs(denom) < 1.0e-30) break;
    const double alpha = rsold / denom;
    for (size_t k = 0; k < n; ++k) {
      x[k] += alpha * pvec[k];
      r[k] -= alpha * Ax[k];
    }
    const double rsnew = dot(r, r);
    if (std::sqrt(rsnew) <= tol) { rsold = rsnew; break; }
    const double beta = rsnew / rsold;
    for (size_t k = 0; k < n; ++k) pvec[k] = r[k] + beta * pvec[k];
    rsold = rsnew;
  }
  common_arrays_C::p = psave;
  const double elapsed =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  std::fprintf(stderr,
               "CPHF nocc=%d nvir=%d |F_ov|rms=%.3e |b|=%.3e iters=%d "
               "resid=%.3e |z|=%.3e time=%.0f ms\n",
               nocc, nvir, fov_rms, bnorm, it, std::sqrt(rsold),
               std::sqrt(dot(x, x)), elapsed);
}
