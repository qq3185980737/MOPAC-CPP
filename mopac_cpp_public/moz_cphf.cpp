// ============================================================================
// moz_cphf.cpp — first building blocks of an exact analytic MOZYME gradient
// (CPHF/LMO-response project).
//
// Stage 1 (this file): reconstruct the dense LMO coefficient matrix C from the
// MOZYME sparse storage (cocc/icocc/ncocc/nncf/ncf), rebuild the AO density
// P = 2 C C^T, and compare it with the packed MOZYME density p.  This validates
// the C layout, which the response module needs.
// Also builds W = C F C^T from the packed AO Fock f (needed for the eventual
// Pulay/response term).
//
// All routines are diagnostics; nothing is called unless an environment
// variable asks for it.
// ============================================================================

#include "moz_cphf.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "ijbo.h"

namespace {
int ijbo_safe(int i, int j) { return ijbo(i, j); }
}  // namespace

namespace MOZYME_C {
extern int noccupied;
}

namespace {

// Build dense C (norbs x nocc, row-major) from the MOZYME LMO storage.
void build_dense_c(std::vector<double>& C) {
  using namespace MOZYME_C;
  const int nr = molkst_C::norbs;
  const int nocc = MOZYME_C::noccupied;
  C.assign((size_t)nr * nocc, 0.0);
  for (int i = 1; i <= nocc; ++i) {
    const int base = ncocc[i];          // starting index of LMO i's atom blocks
    const int nat_i = ncf[i];
    int ja = 0;                          // offset of atom block inside LMO i
    for (int jj = 1; jj <= nat_i; ++jj) {
      const int atom = icocc[nncf[i] + jj];
      const int nj = iorbs[atom];
      const int first = common_arrays_C::nfirst[atom];
      for (int mu = 1; mu <= nj; ++mu) {
        const int row = first - 1 + mu - 1;   // 0-based global AO index
        C[(size_t)row * nocc + (i - 1)] = cocc[base + ja + mu];
      }
      ja += nj;
    }
  }
}

// Packed-density block layout, mirroring helecz/dcart:
//   ijbo(i,j) >= 0 -> 0-based start of the 1-based block, vector index k+1..k+n
//   diagonal block: upper triangle stored (mu >= nu), row-major
//   off-diagonal block (i>j): full iorbs_i x iorbs_j rectangle
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
  if (nr <= 0 || nocc <= 0) {
    std::fprintf(stderr, "MOZCPHF norbs/nocc invalid (%d,%d)\n", nr, nocc);
    return;
  }

  std::vector<double> C;
  build_dense_c(C);

  // P = 2 C C^T
  std::vector<double> P((size_t)nr * nr, 0.0);
  for (int mu = 0; mu < nr; ++mu)
    for (int nu = 0; nu < nr; ++nu) {
      double s = 0.0;
      for (int i = 0; i < nocc; ++i) s += C[(size_t)mu * nocc + i] * C[(size_t)nu * nocc + i];
      P[(size_t)mu * nr + nu] = 2.0 * s;
    }

  // compare with the packed MOZYME density
  std::vector<double> p_dense;
  dense_from_packed(common_arrays_C::p, p_dense);
  double maxerr = 0.0, sum1 = 0.0, sum2 = 0.0;
  for (size_t i = 0; i < P.size(); ++i) {
    maxerr = std::max(maxerr, std::fabs(P[i] - p_dense[i]));
    sum1 += P[i];
    sum2 += p_dense[i];
  }
  std::fprintf(stderr, "MOZCPHF norb=%d nocc=%d |C|_max=%.6f\n", nr, nocc,
               [&] {
                 double m = 0.0;
                 for (double v : C) m = std::max(m, std::fabs(v));
                 return m;
               }());
  std::fprintf(stderr, "MOZCPHF P_dense(max,sum)=(%.10f,%.10f)  P_packed=(%.10f,%.10f)\n",
               [&] { double m = 0.0; for (double v : P) m = std::max(m, std::fabs(v)); return m; }(),
               sum1,
               [&] { double m = 0.0; for (double v : p_dense) m = std::max(m, std::fabs(v)); return m; }(),
               sum2);
  std::fprintf(stderr, "MOZCPHF density max|P_dense-P_packed| = %.3e\n", maxerr);

  // W = C F C^T, with F the packed AO Fock
  std::vector<double> F;
  dense_from_packed(packed_fock, F);
  std::vector<double> FC((size_t)nr * nocc, 0.0);   // F * C
  for (int mu = 0; mu < nr; ++mu)
    for (int i = 0; i < nocc; ++i) {
      double s = 0.0;
      for (int nu = 0; nu < nr; ++nu) s += F[(size_t)mu * nr + nu] * C[(size_t)nu * nocc + i];
      FC[(size_t)mu * nocc + i] = s;
    }
  std::vector<double> W((size_t)nr * nr, 0.0);
  for (int mu = 0; mu < nr; ++mu)
    for (int nu = 0; nu < nr; ++nu) {
      double s = 0.0;
      for (int i = 0; i < nocc; ++i) s += C[(size_t)mu * nocc + i] * FC[(size_t)nu * nocc + i];
      W[(size_t)mu * nr + nu] = s;
    }
  double asym = 0.0, tr = 0.0;
  for (int mu = 0; mu < nr; ++mu) {
    tr += W[(size_t)mu * nr + mu];
    for (int nu = 0; nu < nr; ++nu)
      asym = std::max(asym, std::fabs(W[(size_t)mu * nr + nu] - W[(size_t)nu * nr + mu]));
  }
  std::fprintf(stderr, "MOZCPHF Tr(W) = %.10f   max|W-W^T| = %.3e\n", tr, asym);
}
