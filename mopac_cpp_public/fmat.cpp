// ============================================================================
// 【中文导读】fmat.cpp —— 力常数/能量差分公式（fmat.F90）
// ----------------------------------------------------------------------------
// FORCE/THERMO 的数值二阶导数：对每个笛卡尔变量做 4 阶中心差分
// （±δ、±δ/2 四个位移的梯度），按
//   H(i,j) = [8(g(+δ/2)-g(-δ/2)) - (g(+δ)-g(-δ))] * fact/(12δ)
// 组装力常数矩阵；fact = 4.184/fpc_10*1e21（kcal→erg 换算）。
// 本实现保留同样的差分与单位；dipole/evecs 相关项暂置零（频率/ZPE 正确，
// IR 强度后续补）。缺失处见 docs/TRANSLATION_GAPS.md。
// ============================================================================

#include "fmat.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "molkst_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "dipole.h"
#include "symr.h"
#include "sympop.h"
#include "symmetry_C.h"
#include "molsymy.h"

namespace molkst_C {
extern int numat, nvar;
}
extern void compfg(const std::vector<double>&, bool, double&, bool,
                   std::vector<double>&, bool);

void fmat(std::vector<double>& fmatrx, int& nreal, double, double,
          std::vector<std::vector<double>>& deldip, double,
          std::vector<double>& evecs, bool ts) {
    using namespace molkst_C;
    const int nv = nvar;
    std::fprintf(stderr, "NVAR fmat %d\n", nvar);
    if (nv <= 0) return;

    std::vector<double> x(nv + 1, 0.0);
    for (int i = 1; i <= nv; ++i) x[i] = common_arrays_C::xparam[i];

    std::vector<double> gplus2(nv + 1), gminus2(nv + 1);
    std::vector<double> gplus1(nv + 1), gminus1(nv + 1);
    std::vector<double> grad(nv + 1);
    double escf = 0.0;

    const double fact = 4.184 / funcon_C::fpc_10 * 1.0e21;
    const size_t sz = (size_t)nv * (nv + 1) / 2 + 1;
    if (fmatrx.size() < sz) fmatrx.resize(sz, 0.0);
    for (int i = 1; i < (int)fmatrx.size(); ++i) fmatrx[i] = 0.0;

    if (!ts && std::getenv("SYMR_FMAT")) {
        symr();
        std::fprintf(stderr, "CXXSYMRF nsym %d\n", symmetry_C::nsym);
        std::fprintf(stderr, "CXXSYMRF ipo");
        for (int i = 1; i <= numat; ++i)
            for (int j = 1; j <= std::min(symmetry_C::nsym, 4); ++j)
                std::fprintf(stderr, " %d", symmetry_C::ipo[i][j]);
        std::fprintf(stderr, "\nCXXSYMRF r");
        for (int j = 1; j <= std::min(symmetry_C::nsym, 4); ++j)
            for (int k = 1; k <= 9; ++k) std::fprintf(stderr, " %.6f", symmetry_C::r[k][j]);
        std::fprintf(stderr, "\n");
    }
    const bool precis = molkst_C::keywrd.find(" PREC") != std::string::npos ||
                       molkst_C::keywrd.find(" THERMO_PREC") != std::string::npos;
    std::vector<std::vector<double>> dip_plus(nv + 1, std::vector<double>(4, 0.0));
    std::vector<std::vector<double>> dip_minus(nv + 1, std::vector<double>(4, 0.0));
    std::vector<double> deltas(nv + 1, 1.0 / 120.0);
    auto make_coord = [&](const std::vector<double>& xv) {
        std::vector<std::vector<double>> crd(3, std::vector<double>(numat + 1, 0.0));
        for (int iat = 1; iat <= numat; ++iat)
            for (int k = 1; k <= 3; ++k) crd[k - 1][iat] = xv[(iat - 1) * 3 + k];
        return crd;
    };
    int prev_lu = 0;
    for (int i = 1; i <= nv; ++i) {
        const double save = x[i];
        double delta = 1.0 / 120.0;
        if (precis) {
            // Fortran: determine a good step from the gradient at +delta.
            x[i] = save + delta;
            compfg(x, true, escf, true, grad, true);
            double gnorm = 0.0;
            for (int k = 1; k <= nv; ++k) gnorm += grad[k] * grad[k];
            gnorm = std::sqrt(gnorm);
            if (gnorm > 1.0e-12) delta = delta * 10.0 / gnorm;
            delta = std::min(0.05, std::max(0.005, delta));

            x[i] = save + delta;
            compfg(x, true, escf, true, gplus2, true);
            x[i] = save - delta;
            compfg(x, true, escf, true, gminus2, true);
        }
        // 2-point (default) or the +-delta/2 points of the 4-point formula.
        deltas[i] = delta;
        x[i] = save + 0.5 * delta;
        if (std::getenv("FMAT_TAG")) {
            std::fprintf(stderr, "CXXX %d +", i);
            for (int k = 1; k <= (nv < 9 ? nv : 9); ++k) std::fprintf(stderr, " %.14f", x[k]);
            std::fprintf(stderr, "\n");
        }
        if (std::getenv("FMAT_TAG")) std::fprintf(stderr, "CXXFGBEG %d\n", i);
        compfg(x, true, escf, true, gplus1, true);
        {
            auto crd = make_coord(x);
            std::vector<double> dv(4, 0.0);
            dipole(common_arrays_C::p, crd, dv, 0);
            for (int k = 1; k <= 3; ++k) dip_plus[i][k] = dv[k];
        }
        if (std::getenv("FMAT_TAG")) {
            std::fprintf(stderr, "CXXFG %d +", i);
            for (int k = 1; k <= (nv < 9 ? nv : 9); ++k) std::fprintf(stderr, " %.14f", gplus1[k]);
            std::fprintf(stderr, "\n");
        }
        if (std::getenv("FMAT_G") && i == 1) {
            std::fprintf(stderr, "CXXGRAD");
            for (int k = 1; k <= (nv < 9 ? nv : 9); ++k)
                std::fprintf(stderr, " %.14f", gplus1[k]);
            std::fprintf(stderr, "\n");
        }
        x[i] = save - 0.5 * delta;
        compfg(x, true, escf, true, gminus1, true);
        {
            auto crd = make_coord(x);
            std::vector<double> dv(4, 0.0);
            dipole(common_arrays_C::p, crd, dv, 0);
            for (int k = 1; k <= 3; ++k) dip_minus[i][k] = dv[k];
        }
        if (std::getenv("FMAT_TAG")) {
            std::fprintf(stderr, "CXXFG %d -", i);
            for (int k = 1; k <= (nv < 9 ? nv : 9); ++k) std::fprintf(stderr, " %.14f", gminus1[k]);
            std::fprintf(stderr, "\n");
        }
        x[i] = save;

        auto hval = [&](int k) {
            if (precis) {
                return (8.0 * (gplus1[k] - gminus1[k])
                        - (gplus2[k] - gminus2[k])) / delta * fact / 12.0;
            }
            return (gplus1[k] - gminus1[k]) * 0.5 / delta * fact;
        };
        // Fortran fills the upper block ll..lu (rows 1..i of column i) and
        // then the lower part of column i (rows i..nvar) separately.
        const int ll = prev_lu + 1;
        const int lu = ll + i - 1;
        for (int l = 1; l <= i; ++l) {
            int idx = ll + l - 1;
            fmatrx[idx] += hval(l);
        }
        for (int k = i; k <= nv; ++k) {
            int idx = k * (k - 1) / 2 + i;
            fmatrx[idx] += hval(k);
            if (std::getenv("FMAT_KK"))
                std::fprintf(stderr, "CXXFMATKK %d %d %.14f\n", i, k, fmatrx[idx]);
        }
        prev_lu = lu;
    }

    // Dipole derivatives (F90 fmat): deldip(k,i) = (dip(+)-dip(-))*0.5/delta
    for (int i = 1; i <= nv; ++i)
        for (int k = 1; k <= 3; ++k)
            if ((int)deldip.size() > k && (int)deldip[k].size() > i)
                deldip[k][i] = (dip_plus[i][k] - dip_minus[i][k]) * 0.5 / deltas[i];
    for (size_t i = 0; i < evecs.size(); ++i) evecs[i] = 0.0;
    // nreal is the caller's nvib; do not overwrite it.
}
