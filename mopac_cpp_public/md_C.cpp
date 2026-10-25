// ============================================================================
// 【中文导读】md_C.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 md_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// md_C.cpp — C++ translation of MOPAC 2016 "md_C.F90" (molecular dynamics globals).
#include "md_C.h"

namespace md_C {
std::vector<std::vector<double>> md_vel;
std::vector<std::vector<double>> md_force;
std::vector<std::vector<double>> md_vel_half;

int md_step_max = DEF_MDSTEP;
double dt_md = DEF_DTMD;
double T_target = DEF_TTARGET;
bool l_mdzero = false;
bool md_stop = false;
bool md_test = true;
}  // namespace md_C
