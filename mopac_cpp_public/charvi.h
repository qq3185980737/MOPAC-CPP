// ============================================================================
// 【中文导读】charvi.h —— 对称性接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// charvi.h — C++ translation of MOPAC 2016 "charvi.F90".
#pragma once
#include <vector>

// charvi: character of a vibrational mode under operation ioper.
double charvi(const std::vector<std::vector<double>>& vects, int jorb,
              int ioper, const std::vector<std::vector<double>>& r, int nvecs);
