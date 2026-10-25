// ============================================================================
// 【中文导读】fock2.h —— Fock/CI 超向量接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// fock2.h — C++ translation of MOPAC 2016 "fock2.F90".
#pragma once
#include <vector>

void fock1dorbs(std::vector<double>& f, const std::vector<double>& ptot,
                const std::vector<double>& pa,
                const std::vector<std::vector<double>>& w,
                int& kr, int ia, int ib, int ilim);
void fock2(std::vector<double>& f, const std::vector<double>& ptot,
           std::vector<double>& p, const std::vector<double>& w,
           const std::vector<double>& wj, const std::vector<double>& wk,
           int numat, const std::vector<int>& nfirst,
           const std::vector<int>& nlast, int mode);
