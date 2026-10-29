// surface_hopping_C.h — Tully fewest-switches non-adiabatic MD (C++ port).
//
// Implements overlap-based fewest-switches surface hopping on top of the
// MOPAC C.I. gradient / MD machinery. The non-adiabatic coupling is obtained
// from the overlap of the active-state C.I. eigenvectors at two consecutive
// MD geometries (R_t and R_{t+dt}), i.e. the time-derivative coupling
//   sigma_ab = (O_ab - O_ba) / (2*dt),   O_ab = <psi_a(R_{t+dt})|psi_b(R_t)>,
// which avoids modifying deri1 (no full derivative-coupling vectors needed).
//
// Physical notes / simplifications (documented, not silent):
//   * Real adiabatic states (MOPAC C.I. vectors are real) -> real
//     fewest-switches; the diagonal energy-phase rotation is dropped.
//   * The transported amplitude is the overlap projection
//       c_a = O_{a,i0} (normalized),  i0 = occupied state at R_t.
//   * Hop probability (fewest-switches rate, real form):
//       P_{i->j} = max(0, 2*dt*sigma_ij*c_i*c_j / c_i^2).
//   * On a hop the velocities are rescaled uniformly (scalar lambda) to
//     conserve total energy; upward hops with insufficient kinetic energy are
//     frustrated (rejected). Directional rescaling along the non-adiabatic
//     coupling vector would require deri1-based derivative couplings (future
//     upgrade).
//
// Indexing: state vectors are column-major flat buffers (lab rows x n_act
// cols), 0-based storage, matching meci_C::conf. Energies are in eV (same as
// meci_C::eig). lab is the C.I. dimension.
#pragma once

#include <cstdio>
#include <vector>
#include <string>

namespace sh_C {

// SURF keyword active (parsed by meci.cpp).
extern bool surf_active;
// Number of active (coupled) electronic states (HOP=n keyword, default 3).
extern int n_act;
// C.I. dimension at the last meci() call.
extern int lab_hop;
// Occupied (currently propagated) state, index within the active list (0-based).
extern int i_cur;
// Active-state C.I. eigenvectors at the current geometry R_{t+dt};
// column-major flat (lab_hop*n_act), 0-based. Filled by meci() when SURF.
extern std::vector<double> cur_conf;
// Active-state energies at the current geometry (eV), length n_act.
extern std::vector<double> cur_eig;
// Occupied state's C.I. vector (lab_hop), used by meci() to track the root.
extern std::vector<double> cur_state_vec;
// Phase-gauged active vectors at the previous geometry (R_t of this step),
// carried across MD steps for a consistent non-adiabatic coupling phase.
extern std::vector<double> prev_gauged;
// Cumulative number of accepted hops (diagnostic).
extern double n_hop_acc;
extern double n_hop_att;
// Seed for the hop RNG (reproducible).
extern unsigned long hop_seed;

// Called once by meci() to (re)initialize from the keyword string when needed.
void sh_setup(const std::string& keywrd);

// Called by meci() after diagonalization / root selection (SURF active) to
// fill the active-state data. conf_full is the full C.I. eigenvector matrix
// (column-major lab*lab, 0-based storage); eig_full the full eigenvalues (eV,
// 1-based-padded); active_idx the full-lab 1-based indices of the active
// states in energy order; occupied the full-lab 1-based tracked-state index.
void sh_store_full(const std::vector<double>& conf_full,
                   const std::vector<double>& eig_full,
                   const std::vector<int>& active_idx,
                   int occupied);

// Called by md_driver at the end of each MD step to propagate the electronic
// coefficients, decide whether to hop, and rescale velocities on a hop.
//  conf_prev : active-state C.I. vectors at R_t (flat lab_hop*n_act, 0-based).
//  i0        : occupied-state index within the R_t active ordering.
//  eig_prev  : active-state energies at R_t (eV), length n_act.
//  vel       : md_vel[1..3][1..numat], Bohr / atomic-time units (in/out).
//  mass      : ams[nat[ia]]*U_TO_ME (1-based), in electron-mass units.
//  dt_au     : MD timestep in atomic time units.
//  numat     : number of atoms.
//  iw        : output stream (stdout or a file).
// Returns true if a hop was accepted this step.
bool sh_do_step(const std::vector<double>& conf_prev, int i0,
                const std::vector<double>& eig_prev,
                std::vector<std::vector<double>>& vel,
                const std::vector<double>& mass,
                double dt_au, int numat, FILE* iw);

}  // namespace sh_C
