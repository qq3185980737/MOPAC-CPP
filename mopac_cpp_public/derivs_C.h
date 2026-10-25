// ============================================================================
// 【中文导读】derivs_C.h —— 解析导数/积分接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// derivs_C.h — C++ mapping of Fortran module "derivs_C" (subset used).
#pragma once
#include <vector>
namespace derivs_C {
extern std::vector<double> b;        // derivative scratch
extern std::vector<double> ab;       // (alpha,beta) derivative scratch
extern std::vector<double> fb_ci;    // CI Fock (beta) alias
extern std::vector<double> fmat;     // f matrix scratch
extern std::vector<double> hmat;     // h matrix scratch
extern std::vector<double> wmat;     // w matrix scratch
extern std::vector<double> aidref;   // reference aid array
extern std::vector<double> work2;    // 2-D work array, flat
}
