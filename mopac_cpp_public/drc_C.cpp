// ============================================================================
// 【中文导读】drc_C.cpp —— 反应路径/几何导数
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 drc_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// drc_C.cpp — storage for drc_C.
#include "drc_C.h"
namespace drc_C {
std::vector<double> vref;
std::vector<double> vref0;
std::vector<double> allxyz;
std::vector<double> allvel;
std::vector<double> xyz3;
std::vector<double> vel3;
std::vector<double> allgeo;
std::vector<double> geo3;
std::vector<double> parref;
double time;
}
