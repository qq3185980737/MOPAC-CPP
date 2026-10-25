// ============================================================================
// 【中文导读】geo_diff.h —— 输入/几何接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// geo_diff.h — C++ translation of readmo.F90 "geo_diff" (MOPAC 2016).
#pragma once

// geo_diff calculates the total distance between geometries "geo" and "geoa":
//   sum = addition of all differences in position (sum of motions of atoms).
//   rms = sum of squares of differences in positions.
// If prt is true, lists the atoms that move a lot.
void geo_diff(double& sum, double& rms, bool prt);
