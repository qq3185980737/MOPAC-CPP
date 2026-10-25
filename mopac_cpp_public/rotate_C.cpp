// ============================================================================
// 【中文导读】rotate_C.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 rotate_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// rotate_C.cpp — C++ translation of MOPAC 2016 "rotate_C.F90".
// Fortran uses equivalence ccore(i,1) <-> css_i and ccore(i,2) <-> css2_i;
// here the scalars and the array are kept in sync by callers (ROTATE).
#include "rotate_C.h"

namespace rotate_C {
// thread_local: rotate() runs concurrently under the dcart OpenMP atom-pair
// loop; these Fortran-equivalence scratch globals must be per-thread.
thread_local double css1 = 0.0, csp1 = 0.0, cpps1 = 0.0, cppp1 = 0.0;
thread_local double css2 = 0.0, csp2 = 0.0, cpps2 = 0.0, cppp2 = 0.0;
thread_local double ccore[4][3] = {};  // [4][2] in Fortran; +1 on column for 1-based
}  // namespace rotate_C
