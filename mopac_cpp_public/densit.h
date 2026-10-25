// ============================================================================
// 【中文导读】densit.h —— 密度矩阵接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// densit.h — C++ translation of MOPAC 2016 "densit.F90".
#pragma once
#include <vector>

// Packed lower-triangular density matrix p from square eigenvector matrix c.
void densit(const std::vector<std::vector<double>>& c, int mdim, int norbs,
            int nocc, double occ, int nfract, double fract,
            std::vector<double>& p, int mode);
