// ============================================================================
// 【中文导读】mtxm.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mtxm.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mtxm.cpp
#include "mtxm.h"
void mtxm(const double* a, int nar, const double* b, int nbr, double* c,
          int ncc) {
    for (int i=1;i<=nar;++i)
        for (int j=1;j<=ncc;++j) {
            double s=0;
            for (int k=1;k<=nbr;++k)
                s+=a[(i-1)*nbr+(k-1)]*b[(j-1)*nbr+(k-1)];
            c[(j-1)*nar+(i-1)]=s;
        }
}
