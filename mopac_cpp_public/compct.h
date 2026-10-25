// ============================================================================
// 【中文导读】compct.h —— 结构/成键/电荷接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// compct.h — C++ translation of MOPAC 2016 "compct.F90".
#pragma once
#include <vector>

// Compress LMO data (ic / c / nc / iws) upward, dropping empty LMOs.
void compct(std::vector<int>& nncnew, std::vector<int>& ncnew,
            std::vector<int>& ncmnew, int itop, std::vector<int>& nc,
            std::vector<int>& ic, std::vector<int>& iws, int n01,
            std::vector<double>& c, int n02, int nmos, int idone,
            int& lb, int& mb, int icref, int inref);
