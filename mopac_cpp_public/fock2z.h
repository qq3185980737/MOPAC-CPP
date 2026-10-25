// ============================================================================
// 【中文导读】fock2z.h —— Fock/CI 超向量接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// fock2z.h — C++ translation of MOPAC 2016 "fock2z.F90".
#pragma once
#include <vector>

void focd2z(int iab, int jba,
            std::vector<double>& fii, std::vector<double>& fjj,
            std::vector<double>& fij,
            const std::vector<double>& pii, const std::vector<double>& pjj,
            const std::vector<double>& pij,
            const std::vector<double>& wj, const std::vector<double>& wk,
            bool diagonal, int& kr);
void fock2z(std::vector<double>& f, std::vector<double>& q,
            std::vector<double>& qe, const std::vector<double>& wj,
            const std::vector<double>& wk,
            std::vector<std::vector<double>>& ptot2, int mode, int ione);
