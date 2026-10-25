// ============================================================================
// 【中文导读】dtran2.h —— 解析导数/积分接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// dtran2.h — C++ translation of MOPAC 2016 "dtran2.F90".
#pragma once
#include <vector>

// Transform d-orbital rotation: r(3,3) rotation -> t(5,5,ioper).
void dtran2(std::vector<std::vector<double>>& r,
            std::vector<std::vector<std::vector<double>>>& t, int ioper);
