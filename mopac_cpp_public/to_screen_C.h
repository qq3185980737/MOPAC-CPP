// ============================================================================
// 【中文导读】to_screen_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// to_screen_C.h — C++ mapping of Fortran module "to_screen_C".
#pragma once
#include <vector>
namespace to_screen_C {
// travel(k) - distance traveled during vibration k (1-based).
extern std::vector<double> travel;
// redmas(k,1) reduced mass; redmas(k,2) effective mass.
extern std::vector<std::vector<double>> redmas;
// rot(1..3) rotational constants (cm^-1); xyzmom moments of inertia.
extern std::vector<double> rot;
extern std::vector<double> xyzmom;
extern std::vector<double> freq;    // vibrational frequencies (force.F90)
extern std::vector<double> cnorml;  // normal-mode eigenvectors, packed (force.F90)
extern std::vector<double> dipt;    // vibrational dipole intensities (force.F90)

extern std::vector<std::vector<double>> dip;
}