// ============================================================================
// 【中文导读】supdot.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 supdot.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// supdot.cpp
#include "supdot.h"
#include <cstdio>
// hoff: offset for the packed lower-triangle h. Fortran element n of h is
// stored at physical n-1 when hoff==0 (0-based vectors) and at physical n
// when hoff==1 (1-based-padding buffers like fmat/hmat in deri1).
void supdot(double* s, const double* h, const double* g, int n, int hoff) {
    // hoff: offset for the packed lower-triangle h. Fortran element n of h is
    // stored at physical n-1 when hoff==0 (0-based vectors) and at physical n
    // when hoff==1 (1-based-padding buffers like fmat/hmat in deri1).
    int k=0;
    for (int i=1;i<=n;++i) {
        double sum=0.0;
        for (int j=1;j<=i;++j) {
            sum += g[j]*h[k+j-1+hoff];
        }
        s[i]=sum;
        k+=i;
    }
    if (n==1) return;
    k=1;
    for (int i=2;i<=n;++i) {
        double gi=g[i];
        for (int j=1;j<=i-1;++j) s[j]+=h[k+j-1+hoff]*gi;
        k+=i;
    }
}
