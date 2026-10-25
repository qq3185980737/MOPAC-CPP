// ============================================================================
// 【中文导读】mxm.cpp —— 数值/线性代数与几何工具
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mxm.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mxm.cpp — C++ translation of MOPAC 2016 "mxm.F90".
// Rectangular matrix product C = A*B, all matrices fully packed, column-major.
// Fortran semantics: c(j,i) = sum_k a(k,i)*b(j,k) with 1-based a(nar,nbr),
// b(nbr,ncc), c(nar,ncc).  Callers pass C++ pointers (0-based storage):
// work2.data(), &dxyz[1] (= address of element 1, i.e. base+8), &gradnt[i].
// The 1-based Fortran index p(m,n) therefore maps to ptr[(m-1)*ld + (n-1)].
#include "mxm.h"
void mxm(const double* a, int nar, const double* b, int nbr, double* c,
         int ncc) {
    for (int i = 1; i <= nar; ++i)
        for (int j = 1; j <= ncc; ++j) {
            double s = 0;
            for (int k = 1; k <= nbr; ++k)
                s += a[(k - 1) * nar + (i - 1)] * b[(j - 1) * nbr + (k - 1)];
            c[(j - 1) * nar + (i - 1)] = s;
        }
}
