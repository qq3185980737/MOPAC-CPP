// ============================================================================
// 【中文导读】H_bonds4.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// H_bonds4.h — C++ translation.
#pragma once
#include <vector>
double H_bonds4(bool l_grad, double* dxyz);
double energy_corr_hh_rep(bool l_grad, double* dxyz);
double energy_corr_h4(bool l_grad, double* grad_h4);
double poly(double r, bool l_grad, double& dpoly);
double cvalence_contribution(int atom_a, int atom_b);
double cvalence_contribution_d(int atom_a, int atom_b);
void prt_hbonds(int D, int H, int A, double energy);
