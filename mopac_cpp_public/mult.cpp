// ============================================================================
// 【中文导读】mult.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mult.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mult.cpp — C++ translation of MOPAC 2016 "mult.F90".
// F90: vecs(j,i) = sum_k c(k,i)*s(j,k)  =>  VECS = S*C (column-major).
// Column-major offsets: c(k,i)->(k-1)+(i-1)*n, s(j,k)->(j-1)+(k-1)*n,
//                       vecs(j,i)->(j-1)+(i-1)*n.
#include "mult.h"
void mult(const double* c, const double* s, double* vecs, int n) {
    for (int i=1;i<=n;++i)
        for (int j=1;j<=n;++j) {
            double sum=0;
            for (int k=1;k<=n;++k)
                sum+=c[(k-1)+(i-1)*n]*s[(j-1)+(k-1)*n];
            vecs[(j-1)+(i-1)*n]=sum;
        }
}
