// ============================================================================
// 【中文导读】esp_support.h —— 偶极/极化/静电接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// esp_support.h — helpers for new_esp (MOPAC 2016).
#pragma once
#include <vector>

// S^{-1/2} helper — direct translation of esp_utilities.F90
// "get_minus_point_five_overlap" (rebuilds the normalised overlap from the
// packed h matrix, diagonalises it, forms S^{-1/2}).  Fills s (1-based
// matrix) with S^{-1/2}.
void get_minus_point_five_overlap(std::vector<std::vector<double>>& s);

// Density-matrix construction from MO coefficients (density_for_GPU, GPU
// variant of densit): p(ij) = 2*C_closed*C_closed^T + fract*C_open*C_open^T.
void density_for_GPU(const double* vecs, double fract, int nclose, int nopen,
                     double two, int ij, int norbs, int mode,
                     std::vector<double>& p, int mode2);
