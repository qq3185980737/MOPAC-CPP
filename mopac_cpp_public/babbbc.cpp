// ============================================================================
// 【中文导读】babbbc.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 babbbc.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// babbbc.cpp — C++ translation of MOPAC 2016 "babbbc.F90".
#include "babbbc.h"

#include "meci_C.h"

static inline int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}

double babbbc(const int* iocca1, const int* ioccb1, const int* ioccb2,
              int nmos, const double* xy) {
    int i = 1;
    for (; i <= nmos; ++i) {
        if (ioccb1[i] == ioccb2[i]) continue;
        break;
    }
    int ij = 0;
    int j = i + 1;
    for (; j <= nmos; ++j) {
        if (ioccb1[j] != ioccb2[j]) break;
        ij += iocca1[j] + ioccb1[j];
    }
    ij += iocca1[j];
    double sum = 0.0;
    for (int k = 1; k <= nmos; ++k)
        sum += (xy[xyidx(i, j, k, k, nmos)] - xy[xyidx(i, k, j, k, nmos)]) *
                   (ioccb1[k] - meci_C::occa[k]) +
               xy[xyidx(i, j, k, k, nmos)] * (iocca1[k] - meci_C::occa[k]);
    if ((ij / 2) * 2 != ij) sum = -sum;
    return sum;
}
