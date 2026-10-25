// ============================================================================
// 【中文导读】polar_C.cpp —— 偶极矩/极化率/静电与介电
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 polar_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// polar_C.cpp — C++ translation of MOPAC 2016 "polar_C.F90".
#include "polar_C.h"

namespace polar_C {
double omega = 0.0;    // frequency-dependent polarizability energy (eV)
double alpavg = 0.0;   // average polarizability at omega (cubic Angstroms)
double dummy = 0.0;
}  // namespace polar_C
