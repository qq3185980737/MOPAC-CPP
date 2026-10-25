// ============================================================================
// 【中文导读】chklew.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// chklew.h — C++ translation of MOPAC 2016 "chklew.F90".
#pragma once
#include <vector>

// add_Lewis_element: record a Lewis-structure element and update bookkeeping.
// element_type is incremented for sigma/lone/pi bond types.
void add_Lewis_element(int atom_i, int atom_j, int charge, int& element_type);

// ring5 / arom / arom2: five-membered-ring and aromaticity helpers.
void ring5(int i, const std::vector<int>& mb, std::vector<int>& ir5);
bool arom(int ii, int jj, const std::vector<int>& mpii);
bool arom2(int ii, int jj, const std::vector<int>& mpii);

// chklew: build the Lewis structure (sigma, lone pairs, pi bonds).
void chklew(std::vector<int>& mb, std::vector<int>& numbon, int& l,
            int large, bool debug);
