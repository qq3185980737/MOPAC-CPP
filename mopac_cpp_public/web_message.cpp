// ============================================================================
// 【中文导读】web_message.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 web_message.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// web_message.cpp — print a pointer to the online manual page.
#include "web_message.h"

#include <cstdio>

// 1-arg form used by getgeo.cpp (defaults to the standard output unit).
void web_message(const char* page) {
    web_message(6, page);
}

void web_message(int iw, const char* page) {
    std::printf("\n%10sFor more information, see: HTTP://OpenMOPAC.net/Manual/%s\n\n",
                "", page ? page : "");
}
