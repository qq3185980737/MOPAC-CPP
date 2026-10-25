// ============================================================================
// 【中文导读】diagi.cpp —— CI 能量中的 xy 关联项求值（diagi）
// ----------------------------------------------------------------------------
// 对占据/活性轨道（ialpha/ibeta 过滤）计算 sum_{i,j} xy(i,i,j,j)/xy(i,j,i,j)
// 类型的双电子关联贡献，用于 CI 能量/梯度；xy 为四维打包列主序
// (nmos,nmos,nmos,nmos)。1-based 轨道下标。
// ============================================================================

// diagi.cpp
#include "diagi.h"
double diagi(const int* ialpha, const int* ibeta, const double* eiga,
             const double* xy, int nmos) {
    double x = 0.0;
    for (int i = 1; i <= nmos; ++i) {
        if (ialpha[i] == 0) continue;
        x += eiga[i];
        for (int j = 1; j <= nmos; ++j) {
            // xy(i,i,j,j) and xy(i,j,i,j) in column-major (nmos,nmos,nmos,nmos)
            size_t ij_jj = (size_t)(j - 1) * nmos * nmos * nmos
                         + (j - 1) * nmos * nmos
                         + (i - 1) * nmos
                         + (i - 1);
            size_t ij_ij = (size_t)(j - 1) * nmos * nmos * nmos
                         + (i - 1) * nmos * nmos
                         + (j - 1) * nmos
                         + (i - 1);
            x += (xy[ij_jj] - xy[ij_ij]) * ialpha[j] * 0.5
                 + xy[ij_jj] * ibeta[j];
        }
    }
    for (int i = 1; i <= nmos; ++i) {
        if (ibeta[i] == 0) continue;
        x += eiga[i];
        for (int j = 1; j <= i - 1; ++j) {
            size_t ij_jj = (size_t)(j - 1) * nmos * nmos * nmos
                         + (j - 1) * nmos * nmos
                         + (i - 1) * nmos
                         + (i - 1);
            size_t ij_ij = (size_t)(j - 1) * nmos * nmos * nmos
                         + (i - 1) * nmos * nmos
                         + (j - 1) * nmos
                         + (i - 1);
            x += (xy[ij_jj] - xy[ij_ij]) * ibeta[j];
        }
    }
    return x;
}
