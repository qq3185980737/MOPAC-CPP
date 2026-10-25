// ============================================================================
// 【中文导读】dipole.h —— 偶极/极化/静电接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// dipole.h — C++ translation of MOPAC 2016 "dipole.F90".
#pragma once
#include <vector>

// Dipole moment (Debye). p = packed density, coord(3,n).
double dipole(const std::vector<double>& p,
              std::vector<std::vector<double>>& coord,
              std::vector<double>& dipvec, int mode);
