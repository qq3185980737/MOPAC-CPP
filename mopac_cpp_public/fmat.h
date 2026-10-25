// ============================================================================
// 【中文导读】fmat.h —— Fock/CI 超向量接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// fmat.h — C++ translation of MOPAC 2016 "fmat.F90".
#pragma once
#include <vector>

// 2nd-order central difference of gradient -> Hessian element.
double fmat_second(double grold, double grad, double delta, double fact);
// 4th-order correction (precis).
double fmat_fourth(double grold, double grad, double g2old, double g2rad,
                   double delta, double fact);
void fmat(std::vector<double>& fmatrx, int& nreal, double tscf, double tder,
          std::vector<double>& deldip, double heat,
          std::vector<double>& evecs, bool ts);
