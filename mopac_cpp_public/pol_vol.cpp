// ============================================================================
// 【中文导读】pol_vol.cpp —— CI/激发态与极化
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 pol_vol.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// pol_vol.cpp — C++ translation of "pol_vol.F90".
// Extracted from the full polar response module so that the finite-field
// static_polarizability path does not drag in the TDHF machinery.
#include "pol_vol.h"

#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace parameters_C;

double pol_vol(double average) {
  double polarizability = average * a0 * a0 * a0;
  if (polvol[99] > 1e-4) {
    polarizability = polarizability * polvol[99] + polvol[100];
    for (int i = 1; i <= numat; ++i) polarizability += polvol[labels[i]];
  }
  return polarizability / (a0 * a0 * a0);
}
