// moz_dense_scf.cpp — variational dense-LMO SCF (internal SCF reimplementation).
//
// Opt-in (MOPAC_MOZ_DENSE_SCF=1).  Iteration:
//   P = 2 C C^T -> packed density -> Fock -> residual R = F C - C (C^T F C)
//   C <- orthonormalise(C - eta R) with backtracking line search on the energy
// until max|R| < tol.  No DIIS, no level shift, no local diagonaliser; C is kept
// in static storage across calls (branch memory across geometry changes).
#include "moz_dense_scf.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <vector>

#include "MOZYME_C.h"
#include "buildf.h"
#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "funcon_C.h"
#include "density_for_MOZYME.h"
#include "helecz.h"
#include "makvec.h"
#include "molkst_C.h"
#include "moz_cphf.h"

namespace {

// modified Gram-Schmidt on the columns of C (norbs x n)
void orthonormalise(std::vector<double>& C, int nr, int n) {
  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < j; ++k) {
      double d = 0.0;
      for (int mu = 0; mu < nr; ++mu)
        d += C[(size_t)mu * n + j] * C[(size_t)mu * n + k];
      for (int mu = 0; mu < nr; ++mu) C[(size_t)mu * n + j] -= d * C[(size_t)mu * n + k];
    }
    double nn = 0.0;
    for (int mu = 0; mu < nr; ++mu) {
      const double v = C[(size_t)mu * n + j];
      nn += v * v;
    }
    nn = std::sqrt(nn);
    if (nn < 1.0e-12) nn = 1.0;
    for (int mu = 0; mu < nr; ++mu) C[(size_t)mu * n + j] /= nn;
  }
}

// P = 2 C C^T
void density_from_c(const std::vector<double>& C, int nr, int n, std::vector<double>& Pd) {
  Pd.assign((size_t)nr * nr, 0.0);
  for (int mu = 0; mu < nr; ++mu)
    for (int nu = 0; nu <= mu; ++nu) {
      double s = 0.0;
      for (int i = 0; i < n; ++i) s += C[(size_t)mu * n + i] * C[(size_t)nu * n + i];
      s *= 2.0;
      Pd[(size_t)mu * nr + nu] = s;
      Pd[(size_t)nu * nr + mu] = s;
    }
}

// full Fock + electronic energy at density Pd; returns escf (kcal/mol)
double fock_and_energy(const std::vector<double>& Pd, std::vector<double>& Fd, double& ee) {
  std::vector<double> packed;
  moz_pack_dense(Pd, packed);
  common_arrays_C::p = packed;
  std::vector<double> f(packed.size(), 0.0);
  buildf(f, MOZYME_C::partf, 0);
  ee = helecz();
  double escf = (ee + molkst_C::enuclr) * funcon_C::fpc_9 + molkst_C::atheat;
  if (cosmo_C::useps) escf += cosmo_C::solv_energy * funcon_C::fpc_9;
  moz_unpack_dense(f, Fd);
  return escf;
}

// Jacobi symmetric eigen-decomposition; A is overwritten, columns of V are the
// eigenvectors (unsorted).
void jacobi_eig(std::vector<double>& A, std::vector<double>& V, int n) {
  V.assign((size_t)n * n, 0.0);
  for (int i = 0; i < n; ++i) V[(size_t)i * n + i] = 1.0;
  for (int sweep = 0; sweep < 60; ++sweep) {
    double off = 0.0;
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j) {
        const double v = A[(size_t)i * n + j];
        off += v * v;
      }
    if (off < 1.0e-22) break;
    for (int p = 0; p < n; ++p)
      for (int q = p + 1; q < n; ++q) {
        const double apq = A[(size_t)p * n + q];
        if (std::fabs(apq) < 1.0e-14) continue;
        const double app = A[(size_t)p * n + p];
        const double aqq = A[(size_t)q * n + q];
        const double theta = 0.5 * std::atan2(2.0 * apq, aqq - app);
        const double c = std::cos(theta), sn = std::sin(theta);
        for (int k = 0; k < n; ++k) {
          const double akp = A[(size_t)k * n + p];
          const double akq = A[(size_t)k * n + q];
          A[(size_t)k * n + p] = c * akp - sn * akq;
          A[(size_t)k * n + q] = sn * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          const double apk = A[(size_t)p * n + k];
          const double aqk = A[(size_t)q * n + k];
          A[(size_t)p * n + k] = c * apk - sn * aqk;
          A[(size_t)q * n + k] = sn * apk + c * aqk;
        }
        for (int k = 0; k < n; ++k) {
          const double vkp = V[(size_t)k * n + p];
          const double vkq = V[(size_t)k * n + q];
          V[(size_t)k * n + p] = c * vkp - sn * vkq;
          V[(size_t)k * n + q] = sn * vkp + c * vkq;
        }
      }
  }
}

}  // namespace

bool moz_dense_scf(double& ee_out) {
  using namespace MOZYME_C;
  const int nr = molkst_C::norbs;
  const int nocc = MOZYME_C::noccupied;
  if (nr <= 0 || nocc <= 0) return false;

  static std::vector<double> C;
  if ((int)C.size() != nr * nocc) {
    makvec();
    if (molkst_C::moperr) return false;
    // build the real MOZYME density from the LMOs, then recover C from it
    density_for_MOZYME(common_arrays_C::p, 0, MOZYME_C::noccupied, MOZYME_C::partp);
    MOZYME_C::partp = common_arrays_C::p;
    std::vector<double> Pd0;
    moz_unpack_dense(common_arrays_C::p, Pd0);
    std::vector<double> V;
    jacobi_eig(Pd0, V, nr);
    // eigenvalues are the diagonal of the rotated Pd0; take the top nocc
    std::vector<std::pair<double, int>> ev(nr);
    for (int i = 0; i < nr; ++i) ev[i] = {Pd0[(size_t)i * nr + i], i};
    std::sort(ev.begin(), ev.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });
    C.assign((size_t)nr * nocc, 0.0);
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    for (int j = 0; j < nocc; ++j) {
      const int col = ev[j].second;
      for (int mu = 0; mu < nr; ++mu)
        C[(size_t)mu * nocc + j] = V[(size_t)mu * nr + col] * inv_sqrt2;
    }
    orthonormalise(C, nr, nocc);
  }

  int maxit = 200;
  if (const char* mi = std::getenv("MOPAC_MOZ_DENSE_MAXIT")) maxit = std::atoi(mi);
  double tol = 1.0e-6;
  if (const char* tl = std::getenv("MOPAC_MOZ_DENSE_TOL")) tol = std::atof(tl);

  const bool dense_dbg = (std::getenv("MOPAC_MOZ_DENSE_DBG") != nullptr);
  std::vector<double> Pd, Fd((size_t)nr * nr, 0.0), FC((size_t)nr * nocc, 0.0),
      R((size_t)nr * nocc, 0.0), Ct((size_t)nr * nocc, 0.0);
  if (dense_dbg) {
    const std::vector<double> p_saved = common_arrays_C::p;
    std::vector<double> Pd_old;
    moz_unpack_dense(p_saved, Pd_old);
    std::vector<double> Fd_old;
    double ee_old = 0.0;
    const double escf_old = fock_and_energy(Pd_old, Fd_old, ee_old);
    std::vector<double> p_ours;
    density_from_c(C, nr, nocc, Pd);
    moz_pack_dense(Pd, p_ours);
    double dmax = 0.0;
    if (p_ours.size() == p_saved.size())
      for (std::size_t i = 1; i < p_ours.size(); ++i)
        dmax = std::max(dmax, std::fabs(p_ours[i] - p_saved[i]));
    std::fprintf(stderr, "DENSE_DBG p_from_C_vs_p maxdiff=%.3e escf_old=%.6f trP=%.3f\n",
                 dmax, escf_old, 0.0);
  }
  double ee = 0.0, rmax = 1.0e30;
  int it = 0;
  double escf = 0.0;
  for (it = 1; it <= maxit; ++it) {
    density_from_c(C, nr, nocc, Pd);
    escf = fock_and_energy(Pd, Fd, ee);
    // FC and C^T F C
    for (int mu = 0; mu < nr; ++mu)
      for (int i = 0; i < nocc; ++i) {
        double s = 0.0;
        for (int nu = 0; nu < nr; ++nu) s += Fd[(size_t)mu * nr + nu] * C[(size_t)nu * nocc + i];
        FC[(size_t)mu * nocc + i] = s;
      }
    std::vector<double> Foo((size_t)nocc * nocc, 0.0);
    for (int i = 0; i < nocc; ++i)
      for (int j = 0; j < nocc; ++j) {
        double s = 0.0;
        for (int mu = 0; mu < nr; ++mu)
          s += C[(size_t)mu * nocc + i] * FC[(size_t)mu * nocc + j];
        Foo[(size_t)i * nocc + j] = s;
      }
    rmax = 0.0;
    for (int mu = 0; mu < nr; ++mu)
      for (int i = 0; i < nocc; ++i) {
        double s = FC[(size_t)mu * nocc + i];
        for (int j = 0; j < nocc; ++j) s -= C[(size_t)mu * nocc + j] * Foo[(size_t)j * nocc + i];
        R[(size_t)mu * nocc + i] = s;
        rmax = std::max(rmax, std::fabs(s));
      }
    if (rmax < tol) break;
    // backtracking gradient step
    double alpha = 0.5;
    bool improved = false;
    for (int ls = 0; ls < 8; ++ls) {
      for (size_t k = 0; k < C.size(); ++k) Ct[k] = C[k] - alpha * R[k];
      orthonormalise(Ct, nr, nocc);
      std::vector<double> Pd_t, Fd_t, R_t;
      double ee_t = 0.0;
      density_from_c(Ct, nr, nocc, Pd_t);
      const double escf_t = fock_and_energy(Pd_t, Fd_t, ee_t);
      if (escf_t < escf) {
        C = Ct;
        improved = true;
        break;
      }
      alpha *= 0.5;
    }
    if (!improved) break;
  }

  // leave the module state consistent with the final C
  density_from_c(C, nr, nocc, Pd);
  escf = fock_and_energy(Pd, Fd, ee);
  ee_out = ee;
  std::fprintf(stderr, "DENSE_SCF iters=%d rmax=%.3e escf=%.6f\\n", it, rmax, escf);
  return true;
}
