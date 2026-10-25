// ============================================================================
// 【中文导读】mecid.cpp —— CI/激发态与极化
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mecid.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mecid.cpp - C++ translation of MOPAC 2016 "mecid.F90" + "diagi.F90".
// MECID: remove CI-active electrons from ground state; GSE = stabilization
// energy; EIGA = 1-electron levels; DIAG = microstate energies.
#include "mecid.h"

#include <vector>

#include "meci_C.h"

using namespace meci_C;

// xy(a,b,c,d) column-major index (F90 xy(nmos,nmos,nmos,nmos)).
static int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}

// diagi.F90: energy of a microstate defined by ialpha/ibeta.
static double diagi(const std::vector<int>& ialpha,
                    const std::vector<int>& ibeta, const double* eiga,
                    const double* xy, int nmos) {
    double x = 0.0;
    for (int i = 1; i <= nmos; ++i) {
        if (ialpha[i] == 0) continue;
        x += eiga[i];
        for (int j = 1; j <= nmos; ++j)
            x += (xy[xyidx(i, i, j, j, nmos)] - xy[xyidx(i, j, i, j, nmos)]) *
                     ialpha[j] * 0.5 +
                 xy[xyidx(i, i, j, j, nmos)] * ibeta[j];
    }
    for (int i = 1; i <= nmos; ++i) {
        if (ibeta[i] == 0) continue;
        x += eiga[i];
        for (int j = 1; j < i; ++j)
            x += (xy[xyidx(i, i, j, j, nmos)] - xy[xyidx(i, j, i, j, nmos)]) *
                 ibeta[j];
    }
    return x;
}

void mecid(const double* eigs, double& gse, double* eiga, double* diag,
           double* xy) {
    gse = 0.0;
    for (int i = 1; i <= nmos; ++i) {
        double x = 0.0;
        for (int j = 1; j <= nmos; ++j) {
            double xyij = xy[xyidx(i, i, j, j, nmos)];
            double xyyj = xy[xyidx(i, j, i, j, nmos)];
            x += (2.0 * xyij - xyyj) * occa[j];
        }
        eiga[i] = eigs[i - 1] - x;
        gse += eiga[i] * occa[i] * 2.0;
        gse += xy[xyidx(i, i, i, i, nmos)] * occa[i] * occa[i];
        for (int j = i + 1; j <= nmos; ++j)
            gse += 2.0 * (2.0 * xy[xyidx(i, i, j, j, nmos)] -
                          xy[xyidx(i, j, i, j, nmos)]) *
                   occa[i] * occa[j];
    }
    for (int i = 1; i <= lab; ++i) {
        std::vector<int> ma(nmos + 1), mb(nmos + 1);
        for (int j = 1; j <= nmos; ++j) {
            ma[j] = microa[j][i];
            mb[j] = microb[j][i];
        }
        diag[i] = diagi(ma, mb, eiga, xy, nmos) - gse;
    }
}
