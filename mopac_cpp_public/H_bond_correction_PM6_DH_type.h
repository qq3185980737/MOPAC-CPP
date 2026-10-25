// ============================================================================
// 【中文导读】H_bond_correction_PM6_DH_type.h —— 氢键修正参数声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// H_bond_correction_PM6_DH_type.h — C++ translation.
#pragma once
double PM6_DH_H_bond_corrections(bool l_grad, bool prt);
void setup_DH_Plus(int nrpairs, int* nrbondsa, int* nrbondsb,
                   bool& l_h_bonds, double* covrad);
