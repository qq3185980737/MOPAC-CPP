// ============================================================================
// 【中文导读】maps_C.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 maps_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

#include "maps_C.h"
namespace maps_C {
double rxn_coord=0, rc_escf=0, ekin=0;
int lpara1=0, latom1=0, lpara2=0, latom2=0;
int lparam=0, latom=0, kloop=0;
std::vector<double> react;
double rxn_coord1=0, rxn_coord2=0;
int ione=0, ijlp=0, ilp=0, jlp=0, jlp1=0;
std::vector<double> surf;
}
