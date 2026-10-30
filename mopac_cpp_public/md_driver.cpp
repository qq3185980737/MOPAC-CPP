// ============================================================================
// 闁靛棙鍔掗懙鎴﹀棘閸パ屽殼閻犲洦鐪归埀顒佸敥d_driver.cpp 闁炽儲鏌￠埀?MOPAC 2016 閻犱緤绱曢悾?鐎规悶鍎遍崣鍨熼垾铏仴
// ----------------------------------------------------------------------------
// 閻庣數鎳撶花?MOPAC 2016 md_driver.F90 闁?C++ 缂佸绮ˇ鏌ユ晬閸喎鐏楅柛姘湰缁喚鈧稒鍔楅埢鍏兼償韫囥儳绀嗛柕?// 缂佹拝绠戦悾楣冩晬濮橆偆绠介柟?Fortran 1-based 缂佷究鍨圭槐鈺呭Υ娴ｇ懓鈪甸柛鏍ф噸缁楀懏绋夋径搴健/闁告帗銇炵€靛本鎯旇箛鎾额伌閻忕偐鍋撻柨娑欑☉閸欐寧鎷呴幘宕囨澖闁绘粓顣︾粭?// 闁稿繑濞婇弫顓犫偓娑欘殔閸ㄥ酣寮ㄩ婵愭綄闁哄倸娲ｅ▎銏ゅ礃閸涱厼姣愰柡浣规緲瀵兘宕㈤悢鍛婄畳闁兼槒椴搁弸鍐ㄢ枖閵娾晛娅為柨娑欑〒瀹歌鲸寰?闁告濮崇紞鍛緞閸曨偄鍤掗悹浣规緲缂嶅秹宕?// docs/TRANSLATION_GAPS.md闁?// ============================================================================

// md_driver.cpp 闁?C++ port of the user's Fortran velocity-Verlet MD driver.
//
// Faithful translation of "subroutine run_md_velocity_verlet" and
// "subroutine gen_md_vel" from C:\Users\15632\Desktop\md.txt, operating on the
// existing C++ globals. Two deliberate fixes vs the Fortran original:
//   1. The box size is a real keyword-controlled variable (default 10.0 A),
//      and the printed box always matches the box actually used (the Fortran
//      hard-coded the print as "[0.0, 40.0]" while using 10 A).
//   2. T_target / dt / box can be set by TEMP=, DT=, BOX= keywords instead of
//      only by editing the global defaults.
// Momentum / angular-momentum removal stay disabled, exactly as in the Fortran
// original (they are commented out there).
#include "md_driver.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "compfg.h"
#include "cosmo_C.h"
#include "funcon_C.h"
#include "maps_C.h"
#include "md_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "reada.h"
#include "surface_hopping_C.h"

using namespace common_arrays_C;   // nat, loc, xparam
using namespace cosmo_C;           // solv_energy
using namespace funcon_C;          // fpc_9 (1 eV in kcal/mol)
using namespace maps_C;            // react (initial-velocity vector from VELO)
using namespace md_C;              // md_vel, md_force, md_vel_half, dt_md, T_target, md_step_max, md_stop
using namespace molkst_C;          // numat, nvar, escf, elect, enuclr
using namespace parameters_C;      // ams

// Phase-1: set by iter() when the SCF had to be restarted after itrmax.
// The MD step is then rejected and retried with a smaller dt (see below).
extern bool md_scf_restart;

namespace {
// Atomic-unit conversion constants (identical to the Fortran reference).
constexpr double BOHR_TO_ANG   = 0.529177210903;
constexpr double ANG_TO_BOHR   = 1.0 / BOHR_TO_ANG;
constexpr double HARTREE_TO_EV = 27.21073973;
constexpr double EV_TO_HARTREE = 1.0 / HARTREE_TO_EV;
constexpr double AU_TIME_FS    = 0.024188843265;   // 1 atomic unit of time in fs
constexpr double KB_AU         = 3.166811563e-6;   // Hartree / K (md_C::KB)
constexpr double U_TO_ME       = 1822.888486;      // 1 atomic mass unit in electron masses
}  // namespace

// ---------------------------------------------------------------------------
// gen_md_vel: build the initial MD velocities.
//  - VELO keyword present + velocity block read: use the user's velocities
//    (Angstroms/fs in maps_C::react), then zero the net momentum and rescale
//    to T_target 闁?both DRC-referenced (drc.F90 L290-329, L378-381). This keeps
//    the user's relative directions while setting the temperature.
//  - otherwise: Maxwell-Boltzmann at T_target (Box-Muller), rescaled to T_target.
// ---------------------------------------------------------------------------
static void gen_md_vel(const std::string& keywrd, int iw) {
  md_vel.assign(4, std::vector<double>(numat + 1, 0.0));
  const bool use_velo  = (keywrd.find("VELO") != std::string::npos);
  const bool have_velo = use_velo &&
                         (static_cast<int>(react.size()) >= 3 * numat + 1);

  // Random seed for the Maxwell-Boltzmann cold-atom sampling (ensemble support):
  // each run with a different SEED= gives a different initial velocity draw.
  unsigned seed = 12345u;
  {
    std::size_t ps = keywrd.find("SEED=");
    if (ps != std::string::npos) seed = static_cast<unsigned>(reada(keywrd, static_cast<int>(ps) + 1) + 0.5);
  }

  bool hot_active = false;
  int nhot = 0;
  if (have_velo) {
    // 1 閼?fs = AU_TIME_FS / BOHR_TO_ANG  Bohr/au_time.
    for (int ia = 1; ia <= numat; ++ia)
      for (int c = 1; c <= 3; ++c)
        md_vel[c][ia] = react[(ia - 1) * 3 + c] * AU_TIME_FS / BOHR_TO_ANG;

    // ---- HOT-SPOT (laser-focus energy deposition) scheme --------------
    // Convention: any atom whose VELO velocity is non-zero is a "hot-spot"
    // atom whose given velocity is preserved (it models the energy the laser
    // deposits at the focal point). Atoms whose VELO velocity is zero are
    // "cold": they are initialized to a Maxwell-Boltzmann distribution at
    // T_target and rescaled to exactly T_target. Net-momentum removal and
    // the global rescale are skipped so the hot-spot kinetic energy survives;
    // the overall system temperature then stays low while a few atoms carry
    // the deposited energy that transfers to the rest during the MD.
    std::vector<char> is_hot(numat + 1, 0);
    for (int ia = 1; ia <= numat; ++ia) {
      double v2 = md_vel[1][ia] * md_vel[1][ia] +
                  md_vel[2][ia] * md_vel[2][ia] +
                  md_vel[3][ia] * md_vel[3][ia];
      if (v2 > 1.0e-24) { is_hot[ia] = 1; ++nhot; }
    }
    hot_active = (nhot > 0) && (nhot < numat);

    if (hot_active) {
      std::mt19937 rng(seed);
      std::uniform_real_distribution<double> uni(0.0, 1.0);
      for (int ia = 1; ia <= numat; ++ia) {
        if (is_hot[ia]) continue;
        double mass = ams[nat[ia]] * U_TO_ME;
        double sd   = std::sqrt(KB_AU * T_target / mass);
        for (int c = 1; c <= 3; ++c) {
          double u1 = std::max(1.0e-12, uni(rng));
          u1 = std::min(1.0 - 1.0e-12, u1);
          double u2 = uni(rng);
          double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * 3.14159265358979323846 * u2);
          md_vel[c][ia] = z * sd;
        }
      }
      int ncold = numat - nhot;
      double Ek = 0.0;
      for (int ia = 1; ia <= numat; ++ia) {
        if (is_hot[ia]) continue;
        double mass = ams[nat[ia]] * U_TO_ME;
        Ek += 0.5 * mass * (md_vel[1][ia] * md_vel[1][ia] +
                            md_vel[2][ia] * md_vel[2][ia] +
                            md_vel[3][ia] * md_vel[3][ia]);
      }
      double T_cold = (ncold > 0) ? 2.0 * Ek / (3.0 * ncold * KB_AU) : 0.0;
      if (T_cold > 0.0) {
        double lambda = std::sqrt(T_target / T_cold);
        for (int ia = 1; ia <= numat; ++ia) {
          if (is_hot[ia]) continue;
          for (int c = 1; c <= 3; ++c) md_vel[c][ia] *= lambda;
        }
      }
      double Ek_hot = 0.0;
      for (int ia = 1; ia <= numat; ++ia) {
        if (!is_hot[ia]) continue;
        double mass = ams[nat[ia]] * U_TO_ME;
        Ek_hot += 0.5 * mass * (md_vel[1][ia] * md_vel[1][ia] +
                                md_vel[2][ia] * md_vel[2][ia] +
                                md_vel[3][ia] * md_vel[3][ia]);
      }
      std::fprintf(stdout, "==== hot/cold atom initialisation ====\n");
      std::fprintf(stdout, " cold atoms = %d (Maxwell bath), T = %.1f K\n", ncold, T_target);
      std::fprintf(stdout, " bath: %d cold atoms, T = %.1f K\n", ncold, T_target);
      if (nhot > 0) {
        double T_hot = 2.0 * Ek_hot / (3.0 * nhot * KB_AU);
        std::fprintf(stdout, " hot atoms = %d, T_hot = %10.4f K, Ek_hot = %10.6f Ha\n", nhot, T_hot, Ek_hot);
      }
    } else {
      // No hot atoms (or all atoms hot): original global momentum removal.
      double mx = 0.0, my = 0.0, mz = 0.0, tot = 0.0;
      for (int ia = 1; ia <= numat; ++ia) {
        double m = ams[nat[ia]] * U_TO_ME;
        tot += m; mx += m * md_vel[1][ia]; my += m * md_vel[2][ia]; mz += m * md_vel[3][ia];
      }
      if (tot > 0.0) { mx /= tot; my /= tot; mz /= tot; }
      for (int ia = 1; ia <= numat; ++ia) {
        md_vel[1][ia] -= mx; md_vel[2][ia] -= my; md_vel[3][ia] -= mz;
      }
    }
  } else {
    if (use_velo)
      std::fprintf(stdout,
                   " WARNING: VELO requested but no velocity block read; using random.\n");
    // Maxwell-Boltzmann velocities at T_target (Box-Muller, seed from SEED=).
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    for (int ia = 1; ia <= numat; ++ia) {
      double mass = ams[nat[ia]] * U_TO_ME;          // electron masses
      double sd   = std::sqrt(KB_AU * T_target / mass);
      for (int c = 1; c <= 3; ++c) {
        double u1 = uni(rng), u2 = uni(rng);
        u1 = std::max(1.0e-12, u1);
        u1 = std::min(1.0 - 1.0e-12, u1);
        double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * 3.14159265358979323846 * u2);
        md_vel[c][ia] = z * sd;
      }
    }
  }

  // Global rescale to T_target. Skipped when hot-spot partitioning is active:
  // the cold atoms were already rescaled and the hot atoms must be preserved.
  if (!hot_active) {
    int f_vib = 3 * numat;
    double Ek = 0.0;
    for (int ia = 1; ia <= numat; ++ia) {
      double mass = ams[nat[ia]] * U_TO_ME;
      Ek += 0.5 * mass * (md_vel[1][ia] * md_vel[1][ia] +
                          md_vel[2][ia] * md_vel[2][ia] +
                          md_vel[3][ia] * md_vel[3][ia]);
    }
    double T_calc = 2.0 * Ek / (f_vib * KB_AU);

    std::fprintf(stdout, "==== velocity rescaling ====\n");
    std::fprintf(stdout, " target temperature T_target = %10.4f\n", T_target);
    std::fprintf(stdout, " current temperature T_calc  = %10.4f\n", T_calc);
    if (T_calc > 0.0) {
      double lambda = std::sqrt(T_target / T_calc);
      for (int ia = 1; ia <= numat; ++ia)
        for (int c = 1; c <= 3; ++c) md_vel[c][ia] *= lambda;
    std::fprintf(stdout, " rescale factor lambda = %10.4f\n", lambda);
      Ek = 0.0;
      for (int ia = 1; ia <= numat; ++ia) {
        double mass = ams[nat[ia]] * U_TO_ME;
        Ek += 0.5 * mass * (md_vel[1][ia] * md_vel[1][ia] +
                            md_vel[2][ia] * md_vel[2][ia] +
                            md_vel[3][ia] * md_vel[3][ia]);
      }
      T_calc = 2.0 * Ek / (f_vib * KB_AU);
    std::fprintf(stdout, " final temperature T_final = %10.4f\n", T_calc);
    } else {
    std::fprintf(stdout, " T_calc = 0, skip rescaling\n");
    }
    if (use_velo)
      std::fprintf(stdout,
                   " (闁告帗绻傞～鎰版焻閻斿嘲顔婇柡澶堝劥閸?VELO 闂侇偆鍠庣€规娊宕稿Δ瀣鐎圭寮剁粩濠氭⒖鐠哄搫娅嶉柛鏂诲姂閸ｆ椽鐛懜闈涚樆闁烩晩鍠楅悥锝呫€掗埡浣割唺缂傚倵鏅滈弬?\n");
  } else {
    // Report the overall system temperature after hot-spot partitioning.
    double Ek = 0.0;
    for (int ia = 1; ia <= numat; ++ia) {
      double mass = ams[nat[ia]] * U_TO_ME;
      Ek += 0.5 * mass * (md_vel[1][ia] * md_vel[1][ia] +
                          md_vel[2][ia] * md_vel[2][ia] +
                          md_vel[3][ia] * md_vel[3][ia]);
    }
    double T_sys = 2.0 * Ek / (3.0 * numat * KB_AU);
    std::fprintf(stdout, "==== temperature check: T = %10.4f K, Ek = %10.6f Ha ====\n",
                 T_sys, Ek);
  }
}

// ---------------------------------------------------------------------------
// run_md_velocity_verlet: velocity-Verlet MD with elastic-wall box.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// MOZYME exact (numerical) gradient.
// The factory MOZYME dcart gradient is an approximate gradient: FD checks show
// dE/dx disagreements up to a sign flip (Phase6_MOZYME验证_报告.md, sec. 15).
// This helper evaluates the exact derivative of the MOZYME energy by central
// differences; used for MOZYME MD when MOPAC_MOZ_NUMGRAD=1.
// grad is in kcal/mol/A (the same convention as compfg/deriv).
// ---------------------------------------------------------------------------
static void moz_numeric_gradient(const std::vector<double>& x0, double& escf,
                                 std::vector<double>& grad, int nvar_n,
                                 double h_ang) {
  std::vector<double> gdum(grad.size(), 0.0);
  std::vector<double> xp = x0;
  std::vector<double> xm = x0;
  double ep = 0.0, em = 0.0;
  for (int i = 1; i <= nvar_n; ++i) {
    xp[i] = x0[i] + h_ang;
    xm[i] = x0[i] - h_ang;
    compfg(xp, true, ep, true, gdum, false);
    compfg(xm, true, em, true, gdum, false);
    grad[i] = (ep - em) / (2.0 * h_ang);
    if (i == 1 && std::getenv("MOPAC_MOZ_EPDBG") != nullptr)
      std::fprintf(stderr, "EPDBG i=1 ep=%.10f em=%.10f h=%.6f\n", ep, em, h_ang);
  }
  compfg(x0, true, escf, true, gdum, false);   // restore state at x0
  if (std::getenv("MOPAC_MOZ_EPDBG") != nullptr)
    std::fprintf(stderr, "EPDBG e0=%.10f\n", escf);
  if (std::getenv("MOPAC_MOZ_NUMGRAD_DUMP") != nullptr) {
    for (int gi = 1; gi <= 12 && gi <= nvar_n; ++gi)
      std::fprintf(stderr, "NGNUM %d %.6f\n", gi, grad[gi]);
  }
}

void run_md_velocity_verlet(const std::string& keywrd, int iw) {
  // ---- parse control keywords ----
  // T_target (K), default DEF_TTARGET
  T_target = DEF_TTARGET;
  {
    std::size_t p = keywrd.find("TEMP=");
    if (p == std::string::npos) p = keywrd.find("TTARGET=");
    if (p != std::string::npos) T_target = reada(keywrd, static_cast<int>(p) + 1);
  }
  // dt_md (fs), default DEF_DTMD
  dt_md = DEF_DTMD;
  {
    std::size_t p = keywrd.find("DT=");
    if (p != std::string::npos) dt_md = reada(keywrd, static_cast<int>(p) + 1);
  }
  // number of MD steps, default DEF_MDSTEP
  md_step_max = DEF_MDSTEP;
  {
    std::size_t p = keywrd.find("MDSTEPS=");
    if (p == std::string::npos) p = keywrd.find("CYCLES=");
    if (p != std::string::npos) md_step_max = static_cast<int>(reada(keywrd, static_cast<int>(p) + 1) + 0.5);
  }
  // box [0, box_max] in Angstrom, default 10.0 (matches the Fortran constant)
  double box_min_a = 0.0;
  double box_max_a = 10.0;
  {
    std::size_t p = keywrd.find("BOX=");
    if (p != std::string::npos) box_max_a = reada(keywrd, static_cast<int>(p) + 1);
  }
  const double box_min_bohr = box_min_a * ANG_TO_BOHR;
  const double box_max_bohr = box_max_a * ANG_TO_BOHR;
  // PBC: MOPAC's NDDO force field is *not* periodic (no Ewald / no minimum-image
  // electrostatics).  The old per-atom nearest-image wrap therefore changed the
  // potential discontinuously whenever an atom crossed the box edge: measured
  // +426 kcal/mol over 200 fs in the 120-atom run (9 wrap events, 20-130 kcal
  // each), while non-wrap steps conserved energy to ~0.1 kcal/mol.
  // Phase-0 fix: with PBC= present the system is propagated as an OPEN cluster
  // (no wrap, no wall reflection) so the total energy is conserved exactly.
  // A physically periodic MD would require a periodic Fock/core-core evaluation
  // (future work); elastic walls remain available by simply omitting PBC.
  const bool use_pbc = (keywrd.find("PBC") != std::string::npos);
  if (use_pbc)
    std::fprintf(stdout,
                 " WARNING: PBC requested but the force field is non-periodic;"
                 " periodic wrapping is DISABLED (open boundary) to conserve energy.\n");

  std::vector<double> x_new = xparam;   // 1-based, size nvar+1
  std::vector<double> grad(nvar + 1, 0.0);
  const double dt_nominal = dt_md / AU_TIME_FS;    // fs -> au_time (nominal step)
  double dt = dt_nominal;                          // halved temporarily on SCF retries
  int n_scf_reject = 0;
  int n_jump_reject = 0;      // Phase-2: electronic-state jump rejections
  int jump_retry = 0;
  bool etot_prev_valid = false;
  double etot_prev = 0.0;
  double e_pm7comb = 0.0;   // Phase-4 PM7-COMB-lite pair correction (eV)
  double t_sim_fs = 0.0;

  // ---- helper to extract the (x,y,z) of one atom from x_new (in Angstrom) ----
  auto atom_pos = [&](std::vector<double>& v, int ia) {
    double x = 0.0, y = 0.0, z = 0.0;
    for (int i = 1; i <= nvar; ++i) {
      if (loc[1][i] != ia) continue;
      if (loc[2][i] == 1) x = v[i];
      else if (loc[2][i] == 2) y = v[i];
      else if (loc[2][i] == 3) z = v[i];
    }
    return std::vector<double>{x, y, z};
  };

  // MOZYME guard: the MOZYME (localised-MO) SCF path currently does not provide
  // forces consistent with its energy (measured 200-fs drift -167..-416 kcal vs
  // -2..-10 kcal for the conventional path, dt-independent => force error, not
  // integration error).  Refuse MD so that silently wrong trajectories cannot be
  // produced; the full MOZYME gradient audit is tracked in
  // Phase6_MOZYME验证_报告.md.
  if (keywrd.find(" MOZ") != std::string::npos &&
      std::getenv("MOPAC_ALLOW_MOZYME_MD") == nullptr &&
      std::getenv("MOPAC_MOZ_NUMGRAD") == nullptr) {
    std::fprintf(stdout,
                 " [MD] ERROR: MOZYME is not supported for molecular dynamics in this "
                 "build.\n"
                 " [MD] Measured drift: -0.4..-0.6 kcal/mol/fs even with the "
                 "consistent energy accounting, tight SCF (SCFCRT=1e-8) and "
                 "MOPAC_JUMP_KCAL=3 (conventional path: ~0.02 kcal/mol/fs), i.e. "
                 "-80..-120 kcal over 200 fs.  The MOZYME per-step energy/force "
                 "inconsistency must be fixed before NVE statistics can be trusted "
                 "(see Phase6_MOZYME验证_报告.md).  The factory MOZYME gradient is "
                 "an approximate gradient (FD check), so NVE MD requires either "
                 "MOPAC_MOZ_NUMGRAD=1 (exact numerical gradient, 6N+1 SCF/step, "
                 "energy conserving) or the MOPAC_ALLOW_MOZYME_MD=1 override for "
                 "non-quantitative exploration.\n");
    md_stop = true;
    return;
  }
  // ---- Phase-4 PM7-COMB-lite: H-H / O-O pair corrections (opt-in) ----------
  // Fitted delta-ML v0 corrections (see Phase4_deltaML_v0_闁硅翰鍎遍幉?md):
  //   dE(R) = E_high_level - E_PM7 for H-H (deg-4, 0.5-3.0 A) and O-O (deg-2,
  //   1.0-2.0 A).  Only pair terms are included; O-H is left out because the
  //   purely bond-additive model is inconsistent by ~2 kcal for OH/H2O.
  //   Outside the fitted range the correction is held constant (no force).
  //   Enabled with MOPAC_PM7COMB=1; the default path is unchanged.
  const bool use_pm7comb = (std::getenv("MOPAC_PM7COMB") != nullptr);
  auto pm7comb_pot = [&](const std::vector<double>& xv,
                         std::vector<std::vector<double>>* fout, double& eout) {
    eout = 0.0;
    if (!use_pm7comb) return;
    std::vector<double> xcopy = xv;
    std::vector<std::vector<double>> pos(3, std::vector<double>(numat + 1, 0.0));
    for (int ia = 1; ia <= numat; ++ia) {
      std::vector<double> p = atom_pos(xcopy, ia);
      for (int c = 0; c < 3; ++c) pos[c][ia] = p[c];
    }
    for (int i = 1; i <= numat; ++i) {
      for (int j = i + 1; j <= numat; ++j) {
        const int zi = nat[i], zj = nat[j];
        const bool hh = (zi == 1 && zj == 1);
        const bool oo = (zi == 8 && zj == 8);
        if (!hh && !oo) continue;
        double d[3];
        for (int c = 0; c < 3; ++c) d[c] = pos[c][i] - pos[c][j];
        const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (r < 1.0e-6) continue;
        double e = 0.0, dedr = 0.0;
        // v1.1: absolute correction dE_abs(R) = dE_shape(R) + c0 (c0 makes the
        // asymptotic correction vanish: HOF_ref - HOF_PM7 at the minimum), with
        // a 5th-order smooth switch to zero outside the fitted range.
        auto switch5 = [](double rr, double rin, double rout) -> double {
          if (rr <= rin) return 1.0;
          if (rr >= rout) return 0.0;
          const double x = (rr - rin) / (rout - rin);
          return 1.0 - x * x * x * (10.0 - 15.0 * x + 6.0 * x * x);
        };
        auto dswitch5 = [](double rr, double rin, double rout) -> double {
          if (rr <= rin || rr >= rout) return 0.0;
          const double x = (rr - rin) / (rout - rin);
          return (-30.0 * x * x + 60.0 * x * x * x - 30.0 * x * x * x * x) /
                 (rout - rin);
        };
        if (hh) {
          const double c[5] = {-11.66369646, 37.39194723, -40.85120811,
                               17.22131766, -2.45747941};
          const double c0 = 31.75 / 23.060529;     // +1.377 eV (H2 HOF error)
          const double rin = 2.0, rout = 3.0;
          const double rr = std::min(rout, std::max(0.5, r));
          const double shape = c[0] + c[1] * rr + c[2] * rr * rr +
                               c[3] * rr * rr * rr + c[4] * rr * rr * rr * rr;
          const double sw = switch5(r, rin, rout);
          e = (shape + c0) * sw;
          if (r >= 0.5 && r <= rout) {
            const double dshape = c[1] + 2.0 * c[2] * rr + 3.0 * c[3] * rr * rr +
                                  4.0 * c[4] * rr * rr * rr;
            dedr = dshape * sw + (shape + c0) * dswitch5(r, rin, rout);
          }
        } else {
          const double c[3] = {12.61190992, -13.61407068, 2.75753641};
          const double c0 = 5.45 / 23.060529;      // +0.236 eV (O2 HOF error)
          const double rin = 1.6, rout = 2.2;
          const double rr = std::min(rout, std::max(1.0, r));
          const double shape = c[0] + c[1] * rr + c[2] * rr * rr;
          const double sw = switch5(r, rin, rout);
          e = (shape + c0) * sw;
          if (r >= 1.0 && r <= rout) {
            const double dshape = c[1] + 2.0 * c[2] * rr;
            dedr = dshape * sw + (shape + c0) * dswitch5(r, rin, rout);
          }
        }
        eout += e;                       // eV
        if (fout != nullptr && dedr != 0.0) {
          // force = -dE/dr (eV/A) converted to Ha/bohr
          const double s = -dedr / r * BOHR_TO_ANG / HARTREE_TO_EV;
          for (int c = 0; c < 3; ++c) {
            (*fout)[c + 1][i] += s * d[c];
            (*fout)[c + 1][j] -= s * d[c];
          }
        }
      }
    }
  };

  // ---- print initial coordinates (in Bohr) ----
    std::fprintf(stdout, " recomputing forces at new coordinates\n");
  std::fprintf(stdout, " Atom       X(au)            Y(au)            Z(au)\n");
  for (int ia = 1; ia <= numat; ++ia) {
    std::vector<double> p = atom_pos(x_new, ia);
    std::fprintf(stdout, "%4d %16.6E %16.6E %16.6E\n", ia,
                 p[0] * ANG_TO_BOHR, p[1] * ANG_TO_BOHR, p[2] * ANG_TO_BOHR);
  }
  std::fprintf(stdout, "===============================================\n\n");

  gen_md_vel(keywrd, iw);
  std::fprintf(stdout, ">>>>> END INIT md_vel <<<<<\n");
  std::fprintf(stdout, "\n==================== MD SIMULATION START ====================\n");
  std::fprintf(stdout, " Total MD steps      = %6d\n", md_step_max);
  std::fprintf(stdout, " Timestep DT         = %8.4f au_times\n", dt);
  std::fprintf(stdout, " Target temperature  = %8.2f K\n", T_target);
  std::fprintf(stdout, " Simulation Box      = [%.1f, %.1f] (%s)\n",
               box_min_a, box_max_a,
               use_pbc ? "OPEN (PBC disabled: non-periodic force field)"
                       : "Elastic Wall");
  std::fprintf(stdout, "============================================================\n");

  // ================= Yoshida 4th-order symplectic integrator ==============
  // Opt-in (MOPAC_MD_ORDER4=1): a separate compact path used to verify the
  // dt^4 scaling of the energy drift.  3 force evaluations per MD step.
  // Elastic-wall box only; no PBC, no SURF, no PM7COMB, no jump retry.
  if (std::getenv("MOPAC_MD_ORDER4") != nullptr) {
    const bool o4dbg = (std::getenv("MOPAC_MD_O4DBG") != nullptr);
    const double cbrt2 = std::cbrt(2.0);
    const double w1 = 1.0 / (2.0 - cbrt2);
    const double w0 = -cbrt2 / (2.0 - cbrt2);
    const double wsub[3] = {w1, w0, w1};
    const double f_to_ha = BOHR_TO_ANG / (HARTREE_TO_EV * fpc_9);
    if (std::getenv("MOPAC_PM7COMB") != nullptr)
      std::fprintf(stdout, " [MD] WARNING: PM7COMB is ignored by the 4th-order path.\n");
    std::fprintf(stdout, "\n============ MD (Yoshida 4th-order symplectic) ============\n");
    std::fprintf(stdout, " Total MD steps = %6d, nominal dt = %.4f fs, box = [%.1f, %.1f] A\n",
                 md_step_max, dt_nominal * AU_TIME_FS, box_min_a, box_max_a);
    std::vector<double> grad4(nvar + 1, 0.0);
    auto eval_force4 = [&](const std::vector<double>& xs, double& e_out,
                           std::vector<double>& g_out) {
      if (std::getenv("MOPAC_MOZ_NUMGRAD") != nullptr)
        moz_numeric_gradient(xs, e_out, g_out, nvar, 0.001);
      else
        compfg(xs, true, e_out, true, g_out, true);
    };
    int moz_path = 1;
    if (const char* pth = std::getenv("MOPAC_MOZ_PATH")) moz_path = std::atoi(pth);
    auto eval_force_path = [&](const std::vector<double>& x_from,
                               const std::vector<double>& x_to, double& e_out,
                               std::vector<double>& g_out) {
      if (moz_path > 1) {
        std::vector<double> xi = x_from;
        std::vector<double> gdum(g_out.size(), 0.0);
        double edum = 0.0;
        for (int st = 1; st < moz_path; ++st) {
          const double frac = (double)st / (double)moz_path;
          for (std::size_t i = 1; i < xi.size(); ++i)
            xi[i] = x_from[i] + frac * (x_to[i] - x_from[i]);
          compfg(xi, true, edum, true, gdum, false);
        }
      }
      eval_force4(x_to, e_out, g_out);
    };
    std::vector<std::vector<double>> f4(4, std::vector<double>(numat + 1, 0.0));
    eval_force4(x_new, escf, grad4);
    auto force_from_grad = [&]() {
      for (int i = 1; i <= nvar; ++i)
        f4[loc[2][i]][loc[1][i]] = -grad4[i] * f_to_ha;
    };
    double e_tsym_start = 0.0, e_tsym_turn = 0.0;
    force_from_grad();
    double ek0 = 0.0;
    for (int ia = 1; ia <= numat; ++ia) {
      const double m = ams[nat[ia]] * U_TO_ME;
      for (int c = 1; c <= 3; ++c) ek0 += 0.5 * m * md_vel[c][ia] * md_vel[c][ia];
    }
    e_tsym_start = escf / (fpc_9 * HARTREE_TO_EV) + ek0;
    const double f_vib4 = 3.0 * numat;
    const bool tsym = (std::getenv("MOPAC_O4_TSYM") != nullptr);
    const int nsteps_total = tsym ? 2 * md_step_max : md_step_max;
    for (int istep = 1; istep <= nsteps_total; ++istep) {
      if (tsym && istep == md_step_max + 1) {
        // time reversal point: negate all velocities
        double ek = 0.0;
        for (int ia = 1; ia <= numat; ++ia) {
          const double m = ams[nat[ia]] * U_TO_ME;
          for (int c = 1; c <= 3; ++c) {
            md_vel[c][ia] = -md_vel[c][ia];
            ek += 0.5 * m * md_vel[c][ia] * md_vel[c][ia];
          }
        }
        e_tsym_turn = escf / (fpc_9 * HARTREE_TO_EV) + ek;
      }
      const double e_pot_start = escf / (fpc_9 * HARTREE_TO_EV);
      double ek_start = 0.0;
      for (int ia = 1; ia <= numat; ++ia) {
        const double m0 = ams[nat[ia]] * U_TO_ME;
        for (int c = 1; c <= 3; ++c) ek_start += 0.5 * m0 * md_vel[c][ia] * md_vel[c][ia];
      }
      double work_step = 0.0;
      for (int si = 0; si < 3; ++si) {
        const double dts = dt_nominal * wsub[si];
        for (int ia = 1; ia <= numat; ++ia) {
          const double m = ams[nat[ia]] * U_TO_ME;
          for (int c = 1; c <= 3; ++c) md_vel[c][ia] += f4[c][ia] / m * dts * 0.5;
        }
        const std::vector<double> x_before = x_new;
        for (int i = 1; i <= nvar; ++i) {
          const double dx_bohr = md_vel[loc[2][i]][loc[1][i]] * dts;
          if (o4dbg) work_step += f4[loc[2][i]][loc[1][i]] * dx_bohr;
          x_new[i] += dx_bohr * BOHR_TO_ANG;
        }
        if (!use_pbc) {
          for (int i = 1; i <= nvar; ++i) {
            const int c = loc[2][i];
            if (x_new[i] > box_max_a) {
              x_new[i] = 2.0 * box_max_a - x_new[i];
              md_vel[c][loc[1][i]] = -md_vel[c][loc[1][i]];
            } else if (x_new[i] < box_min_a) {
              x_new[i] = 2.0 * box_min_a - x_new[i];
              md_vel[c][loc[1][i]] = -md_vel[c][loc[1][i]];
            }
          }
        }
        eval_force_path(x_before, x_new, escf, grad4);
        force_from_grad();
        for (int ia = 1; ia <= numat; ++ia) {
          const double m = ams[nat[ia]] * U_TO_ME;
          for (int c = 1; c <= 3; ++c) md_vel[c][ia] += f4[c][ia] / m * dts * 0.5;
        }
      }
      double Ek4 = 0.0;
      for (int ia = 1; ia <= numat; ++ia) {
        const double m = ams[nat[ia]] * U_TO_ME;
        for (int c = 1; c <= 3; ++c)
          Ek4 += 0.5 * m * md_vel[c][ia] * md_vel[c][ia];
      }
      const double e_pot_ha = escf / (fpc_9 * HARTREE_TO_EV);
      const double T4 = 2.0 * Ek4 / (f_vib4 * KB_AU);
      std::fprintf(stdout, "\n ---- MD Step %6d ----\n", istep);
    if (std::getenv("MOPAC_MD_HP") != nullptr)
      std::fprintf(stderr, "MDHP %d %.12f %.12f %.12f\n", istep, e_pot_ha, Ek4,
                   e_pot_ha + Ek4);
      std::fprintf(stdout,
                   " Potential Energy (MD mechanical energy, Hartree) = %12.4f\n", e_pot_ha);
      std::fprintf(stdout, " Total Energy (Potential + Kinetic) = %12.4f\n",
                   e_pot_ha + Ek4);
      std::fprintf(stdout, " Instant Temperature = %10.2f K\n", T4);
      std::fprintf(stdout, " Total Kinetic Energy = %12.6f\n", Ek4);
      if (o4dbg) {
        const double dEpot = e_pot_ha - e_pot_start;
        const double dEk = Ek4 - ek_start;
        std::fprintf(stderr, "O4DBG %d %.10f %.10f %.10f\n", istep,
                     dEpot, work_step, dEpot + dEk);
      }
    }
    std::fprintf(stdout, "\n Actual simulated time = %.3f fs (Yoshida 4th order)\n",
                 md_step_max * dt_nominal * AU_TIME_FS);
    if (std::getenv("MOPAC_O4_TSYM") != nullptr) {
      double ek_end = 0.0;
      for (int ia = 1; ia <= numat; ++ia) {
        const double m = ams[nat[ia]] * U_TO_ME;
        for (int c = 1; c <= 3; ++c)
          ek_end += 0.5 * m * md_vel[c][ia] * md_vel[c][ia];
      }
      const double e_end = escf / (fpc_9 * HARTREE_TO_EV) + ek_end;
      std::fprintf(stderr, "TSYM E_start=%.10f E_turn=%.10f E_end=%.10f\n", e_tsym_start, e_tsym_turn, e_end);
    }
    if (std::getenv("MOPAC_O4_DUMPSTATE") != nullptr) {
      for (int ia = 1; ia <= numat; ++ia) {
        double px = 0.0, py = 0.0, pz = 0.0;
        for (int i = 1; i <= nvar; ++i) {
          if (loc[1][i] != ia) continue;
          const int c = loc[2][i];
          if (c == 1) px = x_new[i];
          else if (c == 2) py = x_new[i];
          else if (c == 3) pz = x_new[i];
        }
        std::fprintf(stderr, "O4STATE %d %.12f %.12f %.12f %.12f %.12f %.12f\n",
                     ia, px, py, pz, md_vel[1][ia], md_vel[2][ia], md_vel[3][ia]);
      }
    }
    xparam = x_new;
    return;
  }

  // velocity Verlet main loop
  // Force caching: the force at the end of a step (R_{t+dt}) is the force at the
  // start of the next step (same R), so it is reused instead of recomputing
  // compfg. This halves the SCF evaluations per step (2x speedup).
  std::vector<double> grad_cache(nvar + 1, 0.0);
  bool have_force_cache = false;
  for (int istep = 1; istep <= md_step_max; ++istep) {
    std::vector<std::vector<int>> hit_flag(4, std::vector<int>(numat + 1, 0));
    // Phase-1: keep a copy of the step-start state so the step can be rejected
    // if the SCF at R_{t+dt} had to be restarted (inconsistent energy/force).
    std::vector<double> x_save = x_new;
    std::vector<std::vector<double>> vel_save = md_vel;
    // SURF scratch: active-state data at R_t, captured after the first compfg.
    std::vector<double> conf_prev;
    std::vector<double> eig_prev;
    int i0 = -1;

    // ---- forces at time t (reuse the cached end-of-previous-step force) ----
    if (have_force_cache && getenv("MOPAC_NOCACHE") == nullptr) {
      grad = grad_cache;
    } else {
      if (std::getenv("MOPAC_MOZ_NUMGRAD") != nullptr) {
        moz_numeric_gradient(x_new, escf, grad, nvar, 0.001);
      } else {
        compfg(x_new, true, escf, true, grad, true);
      }
      if (md_stop) break;
      if (std::getenv("MOPAC_GRADDUMP") != nullptr) {
        static long gradump_count = 0;
        if (gradump_count++ == 0) {
          for (int gi = 1; gi <= nvar; ++gi)
            std::fprintf(stderr, "GRADDMP %d %.12f\n", gi, grad[gi]);
        }
      }
      grad_cache = grad;
    }

    // ---- SURF: capture the R_t active-state data for the coupling ----
    if (sh_C::surf_active && sh_C::lab_hop > 0 &&
        (int)sh_C::cur_eig.size() >= 2 &&
        (int)sh_C::cur_conf.size() ==
            (size_t)sh_C::lab_hop * (int)sh_C::cur_eig.size()) {
      conf_prev = sh_C::cur_conf;
      eig_prev = sh_C::cur_eig;
      i0 = sh_C::i_cur;
      int labh = sh_C::lab_hop;
      // Commit the R_t occupied-state vector so the R_{t+dt} compfg tracks
      // the same electronic state via the overlap.
      if (i0 >= 0 && i0 < (int)sh_C::cur_eig.size() &&
          (int)sh_C::cur_state_vec.size() == labh) {
        for (int q = 0; q < labh; ++q)
          sh_C::cur_state_vec[q] = sh_C::cur_conf[(size_t)i0 * labh + q];
      }
    }

    md_force.assign(4, std::vector<double>(numat + 1, 0.0));
    for (int i = 1; i <= nvar; ++i)
      md_force[loc[2][i]][loc[1][i]] =
          -grad[i] * BOHR_TO_ANG / (HARTREE_TO_EV * fpc_9);
    if (use_pm7comb) {
      double e_dummy = 0.0;
      pm7comb_pot(x_new, &md_force, e_dummy);
    }

    // ---- positions in Bohr ----
    std::vector<std::vector<double>> r_bohr(4, std::vector<double>(numat + 1, 0.0));
    for (int ia = 1; ia <= numat; ++ia) {
      std::vector<double> p = atom_pos(x_new, ia);
      r_bohr[1][ia] = p[0] * ANG_TO_BOHR;
      r_bohr[2][ia] = p[1] * ANG_TO_BOHR;
      r_bohr[3][ia] = p[2] * ANG_TO_BOHR;
    }

    // ---- 1. half-step velocity ----
    md_vel_half.assign(4, std::vector<double>(numat + 1, 0.0));
    for (int ia = 1; ia <= numat; ++ia) {
      double mass = ams[nat[ia]] * U_TO_ME;
      for (int c = 1; c <= 3; ++c)
        md_vel_half[c][ia] = md_vel[c][ia] + md_force[c][ia] / mass * dt * 0.5;
    }

    // ---- 2. position update in Bohr + boundary (elastic-wall reflection or PBC wrap) ----
    for (int ia = 1; ia <= numat; ++ia) {
      double xb = r_bohr[1][ia] + md_vel_half[1][ia] * dt;
      double yb = r_bohr[2][ia] + md_vel_half[2][ia] * dt;
      double zb = r_bohr[3][ia] + md_vel_half[3][ia] * dt;
      if (use_pbc) {
        // Phase-0: open boundary, no wrapping (see warning at setup).
      } else {
        if (xb > box_max_bohr) { xb = box_max_bohr - (xb - box_max_bohr); hit_flag[1][ia] = 1; }
        else if (xb < box_min_bohr) { xb = box_min_bohr + (box_min_bohr - xb); hit_flag[1][ia] = 1; }
        if (yb > box_max_bohr) { yb = box_max_bohr - (yb - box_max_bohr); hit_flag[2][ia] = 1; }
        else if (yb < box_min_bohr) { yb = box_min_bohr + (box_min_bohr - yb); hit_flag[2][ia] = 1; }
        if (zb > box_max_bohr) { zb = box_max_bohr - (zb - box_max_bohr); hit_flag[3][ia] = 1; }
        else if (zb < box_min_bohr) { zb = box_min_bohr + (box_min_bohr - zb); hit_flag[3][ia] = 1; }
      }
      r_bohr[1][ia] = xb; r_bohr[2][ia] = yb; r_bohr[3][ia] = zb;
    }
    // write Bohr coords back to x_new (Angstrom) for the next compfg
    for (int ia = 1; ia <= numat; ++ia) {
      for (int i = 1; i <= nvar; ++i) {
        if (loc[1][i] != ia) continue;
        if (loc[2][i] == 1) x_new[i] = r_bohr[1][ia] * BOHR_TO_ANG;
        else if (loc[2][i] == 2) x_new[i] = r_bohr[2][ia] * BOHR_TO_ANG;
        else if (loc[2][i] == 3) x_new[i] = r_bohr[3][ia] * BOHR_TO_ANG;
      }
    }

    // ---- forces at time t+dt ----
    md_scf_restart = false;
    if (std::getenv("MOPAC_MOZ_NUMGRAD") != nullptr) {
        moz_numeric_gradient(x_new, escf, grad, nvar, 0.001);
      } else {
        compfg(x_new, true, escf, true, grad, true);
      }
    if (md_scf_restart) {
      // The SCF hit itrmax and had to restart; accepting this step injects a
      // spurious 20-130 kcal/mol energy jump.  Reject the step, restore the
      // step-start state and retry the same step with half the timestep.
      x_new = x_save;
      md_vel = vel_save;
      dt *= 0.5;
      --istep;
      ++n_scf_reject;
      std::fprintf(stdout,
                   " [MD] SCF restart at step %d: step rejected, dt -> %.4f fs, retrying\n",
                   istep + 1, dt * AU_TIME_FS);
      if (dt * AU_TIME_FS < 0.01) {
        std::fprintf(stdout,
                     " [MD] dt below 0.01 fs after repeated SCF failures; aborting.\n");
        md_stop = true;
        break;
      }
      continue;
    }
    grad_cache = grad;   // cache for the next step's start force
    have_force_cache = true;   // FIX: cache was never armed -> 2 compfg/step

    std::fprintf(stdout, " full velocity update\n");
    std::fprintf(stdout, " FUNCTION (Heat of formation, kcal/mol) = %16.8E\n", escf);
    std::fprintf(stdout, " Atom       X(A)            Y(A)            Z(A)\n");
    for (int ia = 1; ia <= numat; ++ia) {
      std::vector<double> p = atom_pos(x_new, ia);
      std::fprintf(stdout, "%4d %20.8E %20.8E %20.8E\n", ia, p[0], p[1], p[2]);
    }

    std::vector<std::vector<double>> f_atom(4, std::vector<double>(numat + 1, 0.0));
    for (int i = 1; i <= nvar; ++i)
      f_atom[loc[2][i]][loc[1][i]] =
          -grad[i] * BOHR_TO_ANG / (HARTREE_TO_EV * fpc_9);
    if (use_pm7comb) pm7comb_pot(x_new, &f_atom, e_pm7comb);

    // ---- 3. full velocity update ----
    for (int ia = 1; ia <= numat; ++ia) {
      double mass = ams[nat[ia]] * U_TO_ME;
      for (int c = 1; c <= 3; ++c)
        md_vel[c][ia] = md_vel_half[c][ia] + f_atom[c][ia] / mass * dt * 0.5;
    }
    // elastic wall (only when not PBC): reflect the (already-updated) velocity component
    if (!use_pbc)
      for (int ia = 1; ia <= numat; ++ia)
        for (int c = 1; c <= 3; ++c)
          if (hit_flag[c][ia]) md_vel[c][ia] = -md_vel[c][ia];

    // ---- SURF: electronic propagation + fewest-switches hop decision ----
    if (sh_C::surf_active && i0 >= 0 &&
        (int)conf_prev.size() ==
            (size_t)sh_C::lab_hop * (int)sh_C::cur_eig.size() &&
        (int)eig_prev.size() == (int)sh_C::cur_eig.size()) {
      std::vector<double> mass(numat + 1, 0.0);
      for (int ia = 1; ia <= numat; ++ia) mass[ia] = ams[nat[ia]] * U_TO_ME;
      sh_C::sh_do_step(conf_prev, i0, eig_prev, md_vel, mass, dt, numat, stdout);
    }

    // ---- instantaneous temperature ----
    int f_vib = 3 * numat;
    double Ek = 0.0;
    for (int ia = 1; ia <= numat; ++ia) {
      double mass = ams[nat[ia]] * U_TO_ME;
      Ek += 0.5 * mass * (md_vel[1][ia] * md_vel[1][ia] +
                          md_vel[2][ia] * md_vel[2][ia] +
                          md_vel[3][ia] * md_vel[3][ia]);
    }
    double T_inst = 2.0 * Ek / (f_vib * KB_AU);
    // Phase-6: MD potential-energy accounting.
    // The forces (deriv) include the PM7 post-SCF (dispersion/H-bond) gradients,
    // so the consistent potential is the full compfg energy escf.  For backward
    // compatibility with the validated Phase 2-4 datasets the default keeps the
    // historical (elect+enuclr) accounting; MOPAC_MD_ESCF=1 selects the
    // consistent escf accounting (see Phase6_MOZYME验证_报告.md).
    // Default: the forces include the PM7 post-SCF (dispersion/H-bond) gradients,
    // so the consistent potential energy is the full compfg energy escf.
    // MOPAC_MD_EPOT_OLD=1 restores the historical (elect+enuclr) accounting so
    // that the pre-fix Phase 2-4 datasets can still be reproduced.
    double e_pot_ha = escf / (fpc_9 * HARTREE_TO_EV) + e_pm7comb * EV_TO_HARTREE;
    if (std::getenv("MOPAC_MD_EPOT_OLD") != nullptr)
      e_pot_ha = (elect + enuclr + solv_energy + e_pm7comb) * EV_TO_HARTREE;
    if (std::getenv("MOPAC_MDDBG") != nullptr) {
      // power = sum_i F_i . v_i   [Ha / au_time]
      double power = 0.0;
      for (int ia = 1; ia <= numat; ++ia)
        for (int c = 1; c <= 3; ++c) power += md_force[c][ia] * md_vel[c][ia];
      const double HA_PER_AUTIME_TO_KCAL_PER_FS = 627.5095 / 0.024188843265;
      double power_end = 0.0;
      for (int ia = 1; ia <= numat; ++ia)
        for (int c = 1; c <= 3; ++c) power_end += f_atom[c][ia] * md_vel[c][ia];
      std::fprintf(stderr, "MDDBG %d %.10f %.10f %.10f %.10f %.10f\n", istep,
                   e_pot_ha,
                   power * HA_PER_AUTIME_TO_KCAL_PER_FS,
                   power_end * HA_PER_AUTIME_TO_KCAL_PER_FS, escf,
                   e_pot_ha);
    }

    // ---- energy-jump guard (Phase 2) -------------------------------------
    // A single-determinant SCF can switch to a different electronic solution
    // when the reacting system changes character: the reported total energy
    // then jumps by 10-30 kcal/mol within one step.  Such steps are rejected
    // and retried with half the timestep (max 3 retries, then accepted with a
    // warning so the trajectory is not lost).
    {
      double etot_now = e_pot_ha + Ek;
      if (etot_prev_valid) {
        double dE_kcal = 627.5095 * (etot_now - etot_prev);
        // Threshold is configurable: the corrected potential (MOPAC_PM7COMB) releases
        // much more energy per event, so a smaller threshold may be needed to keep
        // the compensated energy continuous (MOPAC_JUMP_KCAL, default 10 kcal).
        static const double jump_kcal = []() {
          const char* e = std::getenv("MOPAC_JUMP_KCAL");
          return e ? std::atof(e) : 10.0;
        }();
        if (std::fabs(dE_kcal) > jump_kcal) {
          if (jump_retry < 3 && dt * AU_TIME_FS > 0.01) {
            ++n_jump_reject;
            ++jump_retry;
            x_new = x_save;
            md_vel = vel_save;
            dt *= 0.5;
            --istep;
            std::fprintf(stdout,
                         " [MD] energy jump %+.1f kcal at step %d: step rejected,"
                         " dt -> %.4f fs, retrying (%d/3)\n",
                         dE_kcal, istep + 1, dt * AU_TIME_FS, jump_retry);
            continue;
          }
          // Retries exhausted: accept the step, but compensate the electronic-
          // state jump by rescaling the velocities so that the total energy is
          // continuous (the jump is an SCF-solution artefact, not a physical
          // energy source).  Reported T is the compensated one.
          double dE = etot_now - etot_prev;          // Ha
          double Ek_target = Ek - dE;
          std::fprintf(stdout,
                       " [MD] WARNING: energy jump %+.1f kcal at step %d accepted"
                       " after %d retries; ",
                       dE_kcal, istep + 1, jump_retry);
          if (Ek_target > 1.0e-6 && Ek > 1.0e-12) {
            double lam = std::sqrt(Ek_target / Ek);
            for (int ia = 1; ia <= numat; ++ia)
              for (int c = 1; c <= 3; ++c) md_vel[c][ia] *= lam;
            Ek = Ek_target;
            T_inst = 2.0 * Ek / (f_vib * KB_AU);
            etot_now = e_pot_ha + Ek;
            std::fprintf(stdout, "velocities rescaled (lambda=%.4f), energy compensated.\n", lam);
          } else {
            std::fprintf(stdout, "compensation not possible (Ek too small); accepted as-is.\n");
          }
          jump_retry = 0;
        } else {
          jump_retry = 0;
        }
      }
      etot_prev = etot_now;
      etot_prev_valid = true;
    }

    std::fprintf(stdout, "\n ---- MD Step %6d ----\n", istep);
    if (std::getenv("MOPAC_MD_HP") != nullptr)
      std::fprintf(stderr, "MDHP %d %.12f %.12f %.12f\n", istep, e_pot_ha, Ek,
                   e_pot_ha + Ek);
    std::fprintf(stdout, " Potential Energy (MD mechanical energy, Hartree) = %12.4f\n",
                 e_pot_ha);
    std::fprintf(stdout, " Total Energy (Potential + Kinetic) = %12.4f\n",
                 e_pot_ha + Ek);
    std::fprintf(stdout, " Instant Temperature = %10.2f K\n", T_inst);
    std::fprintf(stdout, " Total Kinetic Energy = %12.6f\n", Ek);

    t_sim_fs += dt * AU_TIME_FS;   // accepted step: accumulate real simulated time
    dt = dt_nominal;               // back to nominal dt for the next step

    xparam = x_new;
  }

  std::fprintf(stdout,
               "\n Actual simulated time = %.3f fs (nominal dt = %.4f fs, SCF-rejected steps = %d,"
               " energy-jump rejections = %d)\n",
               t_sim_fs, dt_nominal * AU_TIME_FS, n_scf_reject, n_jump_reject);
  std::fprintf(stdout, "\n==================== MD SIMULATION FINISHED ====================\n");
}





