// ============================================================================
// 【中文导读】to_screen_C.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 to_screen_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// to_screen_C.cpp — storage.
#include "to_screen_C.h"
namespace to_screen_C {
std::vector<double> travel;
std::vector<std::vector<double>> redmas;
std::vector<double> rot(4, 0.0);
std::vector<double> freq;   // vibrational frequencies (force.F90)
std::vector<double> cnorml; // normal-mode eigenvectors, packed (force.F90)
std::vector<double> dipt;   // vibrational dipole intensities (force.F90)
std::vector<double> xyzmom(4, 0.0);

std::vector<std::vector<double>> dip(4, std::vector<double>(3,0.0));
}