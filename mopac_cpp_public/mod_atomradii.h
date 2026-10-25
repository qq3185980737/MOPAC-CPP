// ============================================================================
// 【中文导读】mod_atomradii.h —— MOZYME 接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// mod_atomradii.h — C++ mapping of Fortran module "mod_atomradii".
#pragma once
#include <vector>
namespace mod_atomradii {
// is_metal(z): true for elements treated as metals (negative VDW radius).
extern std::vector<bool> is_metal;
// atom_radius_covalent(z): covalent radius of element z; -1 marks metals.
extern std::vector<double> atom_radius_covalent;
// radius(z): per-atom covalent radius (1-based, filled by setup_mopac_arrays).
extern std::vector<double> radius;
}
