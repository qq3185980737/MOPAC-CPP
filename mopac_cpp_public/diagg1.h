// ============================================================================
// 【中文导读】diagg1.h —— 对角化/本征值接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// diagg1.h — C++ translation of MOPAC 2016 "diagg1.F90".
#pragma once
#include <vector>

// Build significant occupied-virtual LMO matrix elements fmo/ifmo (skeleton).
void diagg1(const std::vector<double>& fao, int nocc, int nvir,
            std::vector<double>& eigv, std::vector<double>& ws,
            std::vector<char>& latoms,
            std::vector<std::vector<int>>& ifmo, std::vector<double>& fmo,
            int fmo_dim, int& nij, int idiagg,
            std::vector<double>& avir, std::vector<double>& aocc,
            std::vector<double>& aov);
