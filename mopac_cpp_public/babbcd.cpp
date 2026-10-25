// ============================================================================
// 【中文导读】babbcd.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 babbcd.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// babbcd.cpp — C++ translation of MOPAC 2016 "babbcd.F90".
#include "babbcd.h"

static inline int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}

double babbcd(const int* iocca1, const int* ioccb1, const int* iocca2,
              const int* ioccb2, int nmos, const double* xy) {
    int ij = 0;
    int i = 1;
    for (; i <= nmos; ++i) {
        if (ioccb1[i] >= ioccb2[i]) continue;
        break;
    }
    int j = i + 1;
    for (; j <= nmos; ++j) {
        if (ioccb1[j] < ioccb2[j]) break;
        ij += iocca2[j] + ioccb2[j];
    }
    ij += iocca2[j];
    int k = 1;
    for (; k <= nmos; ++k) {
        if (ioccb1[k] <= ioccb2[k]) continue;
        break;
    }
    int l = k + 1;
    for (; l <= nmos; ++l) {
        if (ioccb1[l] > ioccb2[l]) break;
        ij += iocca1[l] + ioccb1[l];
    }
    ij += iocca1[l];
    double one = ((ij / 2) * 2 == ij) ? 1.0 : -1.0;
    return (xy[xyidx(i, k, j, l, nmos)] - xy[xyidx(i, l, k, j, nmos)]) * one;
}
