// ============================================================================
// 【中文导读】aabacd.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 aabacd.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// aabacd.cpp — C++ translation of MOPAC 2016 "aabacd.F90" (1-based faithful).
#include "aabacd.h"

static inline int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}

double aabacd(const int* iocca1, const int* ioccb1, const int* iocca2,
              const int* ioccb2, int nmos, const double* xy) {
    int ij = 0, i = 1, j, k, l;
    for (; i <= nmos; ++i) {
        if (iocca1[i] >= iocca2[i]) continue;
        break;
    }
    j = i + 1;
    for (; j <= nmos; ++j) {
        if (iocca1[j] < iocca2[j]) break;
        ij += iocca2[j] + ioccb2[j];
    }
    k = 1;
    for (; k <= nmos; ++k) {
        if (iocca1[k] <= iocca2[k]) continue;
        break;
    }
    l = k + 1;
    for (; l <= nmos; ++l) {
        if (iocca1[l] > iocca2[l]) break;
        ij += iocca1[l] + ioccb1[l];
    }
    ij += ioccb2[i] + ioccb1[k];
    double sum = xy[xyidx(i, k, j, l, nmos)] - xy[xyidx(i, l, k, j, nmos)];
    if (ij % 2 == 1) sum = -sum;
    return sum;
}
