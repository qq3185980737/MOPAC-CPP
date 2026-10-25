// ============================================================================
// 【中文导读】symopr.cpp —— 对称性/特征标
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 symopr.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// symopr.cpp
#include "symopr.h"
void symopr(int numat, double* coord, int jump, double* r) {
    for (int i=1;i<=numat;++i) {
        double help[4] = {0.0, coord[(i-1)*3+0], coord[(i-1)*3+1], coord[(i-1)*3+2]};
        for (int j=1;j<=3;++j) {
            double s=0.0;
            if (jump>=0) {
                for (int k=1;k<=3;++k) s += r[(j-1)*3+(k-1)]*help[k];
            } else {
                for (int k=1;k<=3;++k) s += r[(k-1)*3+(j-1)]*help[k];
            }
            coord[(i-1)*3+(j-1)] = s;
        }
    }
}
