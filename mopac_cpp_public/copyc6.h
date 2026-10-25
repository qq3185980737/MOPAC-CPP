// ============================================================================
// 【中文导读】copyc6.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

#pragma once
#include <vector>
constexpr int N_PARS = 32385*5;
namespace copyc6_data { extern const double pars[N_PARS]; }
// c6ab(iat,jat,iadr,jadr,1..3). maxci(iat) = max compressed index per atom.
void copyc6(int maxc, int max_elem,
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& c6ab,
    std::vector<int>& maxci);
