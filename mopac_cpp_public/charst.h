// ============================================================================
// 【中文导读】charst.h —— 对称性接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// charst.h — C++ translation of MOPAC 2016 "charst.F90".
#pragma once
#include <vector>

// charst: character of CI state istate under operation ioper.
// istate < 0 resets cached microstate state and returns 0.
double charst(const std::vector<std::vector<double>>& vects,
              const std::vector<int>& ntype, int istate, int ioper,
              const std::vector<std::vector<double>>& r, int nvecs, bool& first);
