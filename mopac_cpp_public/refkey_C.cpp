// ============================================================================
// 【中文导读】refkey_C.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 refkey_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// refkey_C.cpp — C++ translation of MOPAC 2016 "refkey_C.F90".
#include "refkey_C.h"

namespace refkey_C {
std::string refkey;  // Fortran: character(len=360), saved
}  // namespace refkey_C
