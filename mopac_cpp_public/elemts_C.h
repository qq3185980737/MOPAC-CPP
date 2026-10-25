// ============================================================================
// 【中文导读】elemts_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// elemts_C.h — C++ mapping of Fortran module "elemts_C".
#pragma once
#include <string>
#include <vector>
namespace elemts_C {
// elemnt(z): chemical symbol (e.g. " C ", " N ") of element z (1-based).
extern std::vector<std::string> elemnt;
// cap_elemnt(z): capitalized element symbol used in PDB labels.
extern std::vector<std::string> cap_elemnt;
// atom_names(z): full atom name per element (len 12).
extern std::vector<std::string> atom_names;
}
