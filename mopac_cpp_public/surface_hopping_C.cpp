// surface_hopping_C.cpp — Tully fewest-switches non-adiabatic MD (C++ port).
//
// See surface_hopping_C.h for the method, conventions and documented
// simplifications. The coupling is obtained from the overlap of the active
// C.I. eigenvectors at two consecutive geometries; no deri1 modification.
#include "surface_hopping_C.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <algorithm>

namespace sh_C {

bool surf_active = false;
int  n_act = 3;
int  lab_hop = 0;
int  i_cur = 0;
std::vector<double> cur_conf;      // lab_hop*n_act, column-major, 0-based
std::vector<double> cur_eig;       // n_act, eV
std::vector<double> cur_state_vec; // lab_hop, occupied-state vector
// Phase-gauged active vectors at the previous geometry (R_t of this step),
// carried across MD steps so the non-adiabatic coupling sigma has a
// consistent phase (the raw adiabatic phases flip sign between geometries).
std::vector<double> prev_gauged;
bool sh_debug = false;
double n_hop_acc = 0.0;
double n_hop_att = 0.0;
unsigned long hop_seed = 1234567UL;

namespace {
constexpr double HARTREE_TO_EV = 27.21073973;
constexpr double EV_TO_HARTREE = 1.0 / HARTREE_TO_EV;

inline double reada(const std::string& line, int istart) {
  // Minimal parser: read the first (signed) number from line[istart..].
  // Mirrors the MOPAC reada semantics closely enough for HOP= keyword parsing.
  int i = istart - 1;
  const int n = (int)line.size();
  while (i < n && (line[i] == ' ' || line[i] == '\t' || line[i] == ',')) ++i;
  int start = i;
  if (i < n && (line[i] == '-' || line[i] == '+')) ++i;
  bool saw = false;
  while (i < n && ((line[i] >= '0' && line[i] <= '9') || line[i] == '.' ||
                   line[i] == 'e' || line[i] == 'E' ||
                   ((line[i] == '-' || line[i] == '+') && saw))) {
    if (line[i] == 'e' || line[i] == 'E') saw = true;
    ++i;
  }
  if (i == start) return 0.0;
  std::string tok = line.substr(start, i - start);
  try { return std::stod(tok); } catch (...) { return 0.0; }
}
}  // namespace

void sh_setup(const std::string& keywrd) {
  surf_active = (keywrd.find("SURF") != std::string::npos);
  sh_debug = (keywrd.find("SHDBG") != std::string::npos);
  n_act = 3;
  std::size_t p = keywrd.find("HOP=");
  if (p != std::string::npos) {
    int v = (int)(reada(keywrd, (int)p + 5) + 0.5);
    if (v >= 1 && v <= 20) n_act = v;
  }
  if (!surf_active) {
    cur_conf.clear();
    cur_eig.clear();
    cur_state_vec.clear();
    prev_gauged.clear();
    lab_hop = 0;
    i_cur = 0;
    n_hop_acc = 0.0;
    n_hop_att = 0.0;
  }
}

// Fill the active-state data from the full C.I. eigensolution. Called by
// meci() when SURF is active, after root selection / diagonalization.
//  conf_full : full eigenvector matrix, column-major lab*lab, 0-based storage.
//  eig_full  : full eigenvalues (eV), length lab+? (1-based-padded).
//  active_idx: full-lab 1-based indices of the active states (energy order).
//  occupied  : full-lab 1-based index of the tracked (occupied) state.
void sh_store_full(const std::vector<double>& conf_full,
                   const std::vector<double>& eig_full,
                   const std::vector<int>& active_idx,
                   int occupied) {
  const int na = (int)active_idx.size();
  const int lab = lab_hop;
  if (na < 1 || lab < 1) return;
  cur_conf.assign((size_t)na * (size_t)lab, 0.0);
  cur_eig.assign((size_t)na, 0.0);
  i_cur = 0;
  for (int k = 0; k < na; ++k) {
    const int s = active_idx[k];
    for (int q = 0; q < lab; ++q)
      cur_conf[(size_t)k * lab + q] = conf_full[(size_t)(s - 1) * lab + q];
    cur_eig[k] = eig_full[s - 1];
    if (s == occupied) i_cur = k;
  }
}

// Overlap-projected transported amplitude and fewest-switches hop decision.
bool sh_do_step(const std::vector<double>& conf_prev, int i0,
                const std::vector<double>& eig_prev,
                std::vector<std::vector<double>>& vel,
                const std::vector<double>& mass,
                double dt_au, int numat, FILE* iw) {
  if (!surf_active) return false;
  const int na = (int)cur_eig.size();
  const int lab = lab_hop;
  if (na < 2 || lab < 1) return false;
  if (i0 < 0 || i0 >= na || i_cur < 0 || i_cur >= na) return false;

  // R_t reference: prefer the phase-gauged vectors carried from the previous
  // step (same geometry as conf_prev, but with a consistent phase). On the
  // first step fall back to the raw conf_prev.
  const std::vector<double>& Rref =
      ((int)prev_gauged.size() == na * lab) ? prev_gauged : conf_prev;

  // Phase-align the current (R_{t+dt}) vectors to the R_t reference: flip each
  // state's sign so its largest overlap with the reference is real-positive.
  // This removes the geometry-to-geometry phase flips that otherwise drive
  // every fewest-switches probability to <=0 (no hops near a crossing).
  std::vector<double> cur_gauged(cur_conf.size(), 0.0);
  for (int a = 0; a < na; ++a) {
    double m = -1.0;
    int bp = 0;
    for (int b = 0; b < na; ++b) {
      double ov = 0.0;
      for (int q = 0; q < lab; ++q)
        ov += cur_conf[(size_t)a * lab + q] * Rref[(size_t)b * lab + q];
      double f = std::fabs(ov);
      if (f > m) { m = f; bp = b; }
    }
    double ovbp = 0.0;
    for (int q = 0; q < lab; ++q)
      ovbp += cur_conf[(size_t)a * lab + q] * Rref[(size_t)bp * lab + q];
    const double sgn = (ovbp >= 0.0) ? 1.0 : -1.0;
    for (int q = 0; q < lab; ++q)
      cur_gauged[(size_t)a * lab + q] = sgn * cur_conf[(size_t)a * lab + q];
  }
  // Remember the gauged current vectors as next step's R_t reference.
  prev_gauged = cur_gauged;

  // ---- 1. overlap matrix O_ab = <psi_a(R_{t+dt}) | psi_b(R_t)> ----
  std::vector<std::vector<double>> O(na, std::vector<double>(na, 0.0));
  for (int a = 0; a < na; ++a)
    for (int b = 0; b < na; ++b) {
      double s = 0.0;
      for (int q = 0; q < lab; ++q)
        s += cur_gauged[(size_t)a * lab + q] * Rref[(size_t)b * lab + q];
      O[a][b] = s;
    }

  // ---- 2. antisymmetric time-derivative coupling sigma_ab ----
  // sigma_ab = (O_ab - O_ba)/(2 dt); note sigma_aa = 0.
  std::vector<std::vector<double>> sig(na, std::vector<double>(na, 0.0));
  for (int a = 0; a < na; ++a)
    for (int b = a + 1; b < na; ++b) {
      double s = (O[a][b] - O[b][a]) / (2.0 * dt_au);
      sig[a][b] = s;
      sig[b][a] = -s;
    }

  // ---- 3. transported amplitude from the R_t occupied state ----
  // c_a = O_{a,i0} (normalized): the amplitude of each R_{t+dt} state that the
  // wavefunction carried on state i0 at R_t projects onto.
  std::vector<double> c(na, 0.0);
  double cn = 0.0;
  for (int a = 0; a < na; ++a) { c[a] = O[a][i0]; cn += c[a] * c[a]; }
  if (cn > 1e-30) {
    double inv = 1.0 / std::sqrt(cn);
    for (int a = 0; a < na; ++a) c[a] *= inv;
  }

  // ---- 4. fewest-switches hopping probability from occupied state ----
  // P_{i->j} = max(0, 2*dt*sigma_ij*c_i*c_j / c_i^2)
  const int i = i_cur;             // occupied at R_{t+dt} (meci tracked it)
  const double ci = c[i];
  const double ci2 = ci * ci + 1e-30;
  std::vector<double> P(na, 0.0);
  for (int j = 0; j < na; ++j) {
    if (j == i) continue;
    double p = 2.0 * dt_au * sig[i][j] * ci * c[j] / ci2;
    P[j] = p > 0.0 ? p : 0.0;
  }

  // DIAG (temporary): per-step coupling/probability/energy trace.
  {
    double maxOffO = 0.0;
    for (int a = 0; a < na; ++a)
      for (int b = 0; b < na; ++b)
        if (a != b) maxOffO = std::max(maxOffO, std::fabs(O[a][b]));
    double maxP = 0.0;
    for (int j = 0; j < na; ++j) maxP = std::max(maxP, P[j]);
    if (iw) {
      if (sh_debug) {
        std::fprintf(iw, " SHDBG i0=%d i=%d maxOffO=%.4f maxP=%.5f c[%d]=%.4f eigs=",
                     i0, i, maxOffO, maxP, i, ci);
        for (int j = 0; j < na; ++j) std::fprintf(iw, "%.4f ", cur_eig[j]);
        std::fprintf(iw, "\n");
      }
    }
  }

  // ---- 5. stochastic hop decision ----
  int hop_to = -1;
  double cum = 0.0;
  static std::mt19937 rng(hop_seed);
  double r = std::uniform_real_distribution<double>(0.0, 1.0)(rng);
  for (int j = 0; j < na; ++j) {
    cum += P[j];
    if (cum > 1.0) cum = 1.0;
    if (r < cum) { hop_to = j; break; }
  }
  if (hop_to < 0 || hop_to == i) return false;

  // ---- 6. energy-conserving velocity rescaling ----
  // dE (eV) = E_j - E_i at R_{t+dt}; convert to Hartree.
  const double dE_Ha = (cur_eig[hop_to] - cur_eig[i]) * EV_TO_HARTREE;
  // Kinetic energy (Hartree) = 1/2 sum m_a |v_a|^2.
  double KE = 0.0;
  for (int ia = 1; ia <= numat; ++ia)
    for (int c = 1; c <= 3; ++c)
      KE += 0.5 * mass[ia] * vel[c][ia] * vel[c][ia];
  double KE_new = KE - dE_Ha;
  if (KE_new < 0.0) {
    // Frustrated hop: not enough kinetic energy to climb onto the upper state.
    n_hop_att += 1.0;
    if (iw)
      std::fprintf(iw, " SH: frustrated hop %d->%d (dE=%.4f eV > KE=%.4f Ha); rejected\n",
                   i + 1, hop_to + 1, cur_eig[hop_to] - cur_eig[i], KE);
    return false;
  }
  const double lambda = (KE > 1e-30) ? std::sqrt(KE_new / KE) : 0.0;
  for (int ia = 1; ia <= numat; ++ia)
    for (int c = 1; c <= 3; ++c)
      vel[c][ia] *= lambda;

  // ---- 7. accept the hop; collapse onto the new state ----
  ++n_hop_att;
  ++n_hop_acc;
  i_cur = hop_to;
  cur_state_vec.assign(lab, 0.0);
  for (int q = 0; q < lab; ++q)
    cur_state_vec[q] = cur_conf[(size_t)hop_to * lab + q];
  if (iw)
    std::fprintf(iw, " SH: hop %d->%d accepted (dE=%.4f eV, lambda=%.4f, KE %.4f->%.4f Ha)\n",
                 i + 1, hop_to + 1, cur_eig[hop_to] - cur_eig[i], lambda,
                 KE, KE_new);
  return true;
}

}  // namespace sh_C
