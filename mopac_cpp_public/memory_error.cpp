// ============================================================================
// 【中文导读】memory_error.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 memory_error.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// memory_error.cpp — C++ translation of the runtime support routine.
#include "memory_error.h"

#include <cstdio>
#include <cstdlib>

void memory_error(const char* where) {
    std::printf(" MOPAC  Memory allocation failure in %s\n", where ? where : "?");
    std::fflush(stdout);
    std::exit(1);
}
