// ============================================================================
// 【中文导读】epseta.cpp —— 偶极矩/极化率/静电与介电
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 epseta.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// epseta.cpp — C++ translation of MOPAC 2016 "epseta.F90".

#include "epseta.h"

#include <algorithm>
#include <cfloat>
#include <cmath>

void epseta(double& eps, double& eta) {
    eps = std::max(DBL_MIN, 1e-39);
    eta = std::max(DBL_EPSILON, 1e-17);
}
