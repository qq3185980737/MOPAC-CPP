// ============================================================================
// 【中文导读】rotmol.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 rotmol.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// rotmol.cpp
#include "rotmol.h"
extern void symopr(int numat, double* coord, int, double* r);
void rotmol(int numat, double* coord, double sina, double cosa, int i,
            int j, double* r) {
    symopr(numat, coord, -1, r);
    for (int k=1;k<=3;++k) {
        double buff = -sina*r[(i-1)*3+(k-1)] + cosa*r[(j-1)*3+(k-1)];
        r[(i-1)*3+(k-1)] = cosa*r[(i-1)*3+(k-1)] + sina*r[(j-1)*3+(k-1)];
        r[(j-1)*3+(k-1)] = buff;
    }
    symopr(numat, coord, 1, r);
}
