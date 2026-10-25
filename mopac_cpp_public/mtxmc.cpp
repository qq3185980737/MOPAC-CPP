// ============================================================================
// 【中文导读】mtxmc.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mtxmc.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mtxmc.cpp
#include "mtxmc.h"
#include "mxm.h"
// C(NAR,NAR) = (A(NBR,NAR))' * B(NBR,NAR), packed lower triangle
void mtxmc(const double* a, int nar, const double* b, int nbr, double* c) {
    int l = 0;
    for (int i = 1; i <= nar; ++i) {
        mxm(a + (i-1)*nbr, 1, b, nbr, c + l, i);
        l += i;
    }
}
