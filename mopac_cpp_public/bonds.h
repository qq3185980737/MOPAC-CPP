// ============================================================================
// 【中文导读】bonds.h —— 结构/成键/电荷接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// bonds.h — C++ translation of MOPAC 2016 "bonds.F90".
#pragma once
#include <vector>

// bonds: compute and print bond orders and valencies.
void bonds();

// dopen: ROHF open-shell density matrix from eigenvector matrix c(mdim,mdim).
// sdm(lower triangle, packed) = frac * sum_{nl1..nu1} c(i,n)*c(j,n).
void dopen(const std::vector<std::vector<double>>& eigvec, int mdim, int nor,
           int ndubl, int nsingl, double frac_in, std::vector<double>& sdm);
