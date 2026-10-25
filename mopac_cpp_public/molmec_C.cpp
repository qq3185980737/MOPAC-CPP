// ============================================================================
// 【中文导读】molmec_C.cpp —— MOZYME 局域轨道/成键
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 molmec_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// molmec_C.cpp — C++ translation of MOPAC 2016 "molmec_C.F90".
#include "molmec_C.h"

namespace molmec_C {
int nhco[4][4001] = {};
int nnhco = 0;
double htype = 0.0;
}  // namespace molmec_C
