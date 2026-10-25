// ============================================================================
// 【中文导读】pol_vol.h —— CI/激发态接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// pol_vol.h — C++ translation of "pol_vol.F90".
#pragma once
// pol_vol: average polarizability in A^3 from the mean eigenvalue (a.u.),
// scaled by the "POLVOL" keyword overrides if present.
double pol_vol(double average);
