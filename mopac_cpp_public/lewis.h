// ============================================================================
// 【中文导读】lewis.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// lewis.h — C++ translation of lewis.F90 (MOPAC 2016).
#pragma once
#include <string>

// Build the Lewis structure: raw topography from set_up_dentate is tidied by
// filters (max bonds, hydrogen-bond checks, nitro-group O-O bridge, CVB).
// On exit nbonds/ibonds hold the final connectivity (max 4 bonds per atom).
void lewis(bool use_cvs);

// Remove the longest bond from atom i (used to cap Fe/Ni/Pd/Pt at 5 bonds).
void remove_bond(int i);

// Check hydrogen bonding: nearest non-H neighbour, remove bridge bonds,
// report bad positions (ibad count).  Only active for MOZYME runs.
void check_h(int& ibad);

// Process the CVB keyword: add/delete bonds defined as CVB(j,l) or
// CVB("label","label"); verifies geometry and atom labels.
void check_CVS(bool let);

// Convert a PDB/Jmol atom label in "text" (1-based position j_in after the
// opening quote) to an atom number.  On success text is rewritten with the
// number.  m receives the compressed label length (or the Jmol ']' position).
void txt_to_atom_no(std::string& text, int j_in, bool let, int& m);
