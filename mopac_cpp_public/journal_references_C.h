// ============================================================================
// 【中文导读】journal_references_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// journal_references_C.h — C++ translation of MOPAC 2016 "journal_references_C.F90".
// Per-element journal reference strings for each Hamiltonian. Fortran uses
// equivalence (refmn..refrm1) <-> allref(107,7); here the 7 arrays are kept
// independent (allref view omitted; semantic content identical).
#pragma once

#include <string>

namespace journal_references_C {
// 1-based [108]; entries without a data statement stay "" (e.g. refmn/refam/
// refpm3 are filled at runtime by the parameter loader in Fortran).
extern std::string refmn[108];
extern std::string refpm6[108];
extern std::string refam[108];
extern std::string refpm3[108];
extern std::string refmd[108];
extern std::string refrm1[108];
extern std::string refpm7[108];
}  // namespace journal_references_C
