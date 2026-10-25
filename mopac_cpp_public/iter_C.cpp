// ============================================================================
// 【中文导读】iter_C.cpp —— SCF 迭代/收敛/Fock 构建
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 iter_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// iter_C.cpp — storage.
#include "iter_C.h"
namespace iter_C {
std::vector<double> pold;
}
#include "iter_C.h"
namespace iter_C { std::vector<double> pold2,pold3,pbold,pbold2,pbold3,vec_ai,fock_ai,p_ai,h_ai,vecl_ai,vec_bi,fock_bi,p_bi,h_bi,vecl_bi,theta; double eold_alpha=0,eold_beta=0; bool md_test=false,is_PARAM=false,lgpu=false; }