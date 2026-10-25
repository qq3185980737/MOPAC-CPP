// ============================================================================
// 【中文导读】schmib.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 schmib.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// schmib.cpp
#include "schmib.h"
#include <cmath>
void schmib(double* u, int n, int ndim) {
    const double small=0.01, one=1.0;
    int n1=n+1, ii=0;
    for (int k=1;k<=n;++k) {
        int k1=k-1, col=n1-k;
        double dot=0.0;
        for (int i=1;i<=n;++i) { double v=u[(i-1)+(col-1)*ndim]; dot+=v*v; }
        if (std::fabs(dot)<1e-20) goto l100;
        { double scale=one/std::sqrt(dot);
          for (int i=1;i<=n;++i) u[(i-1)+(col-1)*ndim]*=scale; }
      l30:
        if (k1==0) continue;
        int npass=0;
      l40:
        ++npass;
        for (int j=1;j<=k1;++j) {
            int cj=n1-j;
            double d=0.0;
            for (int i=1;i<=n;++i) d+=u[(i-1)+(cj-1)*ndim]*u[(i-1)+(col-1)*ndim];
            for (int i=1;i<=n;++i) u[(i-1)+(col-1)*ndim]-=d*u[(i-1)+(cj-1)*ndim];
        }
        dot=0.0;
        for (int i=1;i<=n;++i) { double v=u[(i-1)+(col-1)*ndim]; dot+=v*v; }
        if (std::fabs(dot)<1e-20) goto l100;
        if (dot<small && npass>2) goto l100;
        { double scale=one/std::sqrt(dot);
          for (int i=1;i<=n;++i) u[(i-1)+(col-1)*ndim]*=scale; }
        if (dot<small) goto l40;
        continue;
      l100:
        ++ii;
        u[(ii-1)+(col-1)*ndim]=one;
        goto l30;
    }
}
