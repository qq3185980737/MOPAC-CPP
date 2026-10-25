// ============================================================================
// 【中文导读】atomrs.h —— 结构/成键/电荷接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// atomrs.h — C++ translation of MOPAC 2016 "atomrs.F90".
#pragma once
#include <string>
#include <vector>

// atomrs: identify a protein residue starting at iatom, walk the bond graph,
// determine the residue type, and label atoms via txtatm.
// Arrays keep Fortran 1-based indexing (index 0 is padding).
void atomrs(std::vector<int>& lused, std::vector<bool>& ioptl, int& ires,
            int n1, int io, int uni_res, bool first_res);

// peptide_n(l): true if atom l looks like a backbone peptide nitrogen
// (three bonds: two C, one H, with one C attached to a terminal O).
bool peptide_n(int l);
