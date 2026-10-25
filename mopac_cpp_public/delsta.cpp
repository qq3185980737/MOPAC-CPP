// delsta.cpp — C++ translation of MOPAC 2016 "delsta.F90".
//
// ============================================================================
// 【中文导读】delsta() —— MOZYME 局部电荷/密度差分辅助
// ----------------------------------------------------------------------------
// 对应 delsta.F90。用 ijbo 定位原子对在 H/P/F 中的块起点，计算局部电荷
// 差分；2026-10 已删除局部空 ijbo，改调真实实现。
// ============================================================================

#include "delsta.h"

#include <cmath>

#include "funcon_C.h"
#include "molkst_C.h"
#include "ijbo.h"
#include "parameters_C.h"

using namespace funcon_C;
using namespace molkst_C;
using namespace parameters_C;

void delsta(const std::vector<int>& nat, const std::vector<int>& iorbs,
            const std::vector<double>& p,
            const std::vector<std::vector<double>>& cdi,
            std::vector<double>& dstat, int ii, int jj) {
    double qii = tore[nat[ii]];
    int l = ijbo(ii, ii);
    for (int k = 1; k <= iorbs[ii]; ++k) { l += k; qii -= p[l]; }
    double qjj = tore[nat[jj]];
    l = ijbo(jj, jj);
    for (int k = 1; k <= iorbs[jj]; ++k) { l += k; qjj -= p[l]; }
    double rij = std::sqrt(std::pow(cdi[1][1] - cdi[1][2], 2) +
                           std::pow(cdi[2][1] - cdi[2][2], 2) +
                           std::pow(cdi[3][1] - cdi[3][2], 2));
    std::vector<double> vect(4, 0.0);
    if (rij > cutofp) {
        for (int k = 1; k <= 3; ++k) dstat[k] = 0.0;
    } else {
        for (int k = 1; k <= 3; ++k) vect[k] = (cdi[k][1] - cdi[k][2]) / rij;
        double sum = fpc_9 * ev / (rij * rij);
        for (int k = 1; k <= 3; ++k) dstat[k] = -0.5 * qjj * qii * sum * vect[k];
    }
}
