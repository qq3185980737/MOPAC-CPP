// ============================================================================
// 【中文导读】chkion.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// chkion.h — C++ translation of MOPAC 2016 "chkion.F90".
#pragma once
#include <string>
#include <vector>

// chkion: determine which atoms are ionized from Lewis structure.
// atom_charge(i): user-specified charge override (' ', '+', '-', '0').
void chkion(std::vector<int>& ox_calc, int& n_lone_pairs,
            const std::vector<char>& atom_charge);
