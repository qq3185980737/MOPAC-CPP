// ============================================================================
// 【中文导读】quadr.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 quadr.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// quadr.cpp — solve quadratic through three points.
#include "quadr.h"
void quadr(double f0, double f1, double f2, double x1, double x2,
           double& a, double& b, double& c) {
    c = (x2*(f1-f0) - x1*(f2-f0))/(x2*x1*(x1-x2));
    b = (f1 - f0 - c*x1*x1)/x1;
    a = f0;
}
