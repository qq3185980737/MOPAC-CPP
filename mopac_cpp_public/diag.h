// ============================================================================
// 【中文导读】diag.h —— 对角化/本征值接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// diag.h — C++ translation of MOPAC 2016 "diag.F90".
#pragma once
#include <vector>

// Fast pseudo-diagonalization: rotate old eigenvectors to block-diagonalize.
void diag(const std::vector<double>& fao,
          std::vector<std::vector<double>>& vector,
          int nocc, const std::vector<double>& eig, int mdim, int n);
