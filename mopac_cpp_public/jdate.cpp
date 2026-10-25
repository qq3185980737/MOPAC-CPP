// ============================================================================
// 【中文导读】jdate.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 jdate.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// jdate.cpp
#define _CRT_SECURE_NO_WARNINGS
#include "jdate.h"
#include <ctime>
#include <cstdio>

std::string jdate() {
    static const int dim[12] = {0,31,28,31,30,31,30,31,31,30,31,30};
    std::time_t t = std::time(nullptr);
    std::tm* tm = std::localtime(&t);
    int iyear = tm->tm_year % 100;
    int imonth = tm->tm_mon + 1;
    int iday = tm->tm_mday;
    int ijulian = iday;
    for (int i = 1; i <= imonth; ++i) ijulian += dim[i-1];
    char buf[9];
    std::snprintf(buf, sizeof(buf), "%02d%03d", iyear, ijulian);
    return std::string(buf) + "   ";
}
