// ============================================================================
// 【中文导读】H_bond_correction_EH_plus.h —— 氢键修正参数声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// H_bond_correction_EH_plus.h — C++ translation.
#pragma once
double EH_plus(int i, const int* hblist, int max_h_bonds,
               const int* nrbondsa, const int* nrbondsb);
