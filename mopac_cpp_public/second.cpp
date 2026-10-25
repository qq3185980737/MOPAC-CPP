// ============================================================================
// 【中文导读】second.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 second.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// second.cpp
#include "second.h"
#include <chrono>
static std::chrono::steady_clock::time_point t0;
static bool initialized=false;
double second(int) {
    if (!initialized) { t0=std::chrono::steady_clock::now(); initialized=true; }
    auto t1=std::chrono::steady_clock::now();
    return std::chrono::duration<double>(t1-t0).count();
}
