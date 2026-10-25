// ============================================================================
// 【中文导读】timer.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 timer.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// timer.cpp
#include "timer.h"
#include <cstdio>
double CPU_0 = 0.0, CPU_1 = 0.0;
double wall_clock_0 = 0.0, wall_clock_1 = 0.0;
extern double second(int);
void timer(const std::string& a) {
    static bool first = true;
    static double CPU_start, CPU_last;
    static double wall_clock_last;
    if (first) {
        CPU_start = CPU_0;
        CPU_last = CPU_start;
        wall_clock_last = wall_clock_1;
        wall_clock_0 = wall_clock_last;
        first = false;
    }
    int arg1 = 1;
    double dummy = second(arg1);
    if (dummy < -300.0) return;
    double CPU_now = CPU_1;
    std::printf("%s CPU INTERVAL: %.2f, INTEGRAL: %.2f WALL-CLOCK INTERVAL: %.2f, INTEGRAL: %.2f\n",
        a.c_str(), CPU_now - CPU_last, CPU_now - CPU_start,
        wall_clock_1 - wall_clock_last, wall_clock_1 - wall_clock_0);
    CPU_last = CPU_now;
    wall_clock_last = wall_clock_1;
}
