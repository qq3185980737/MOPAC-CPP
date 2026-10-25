// ============================================================================
// 【中文导读】overlaps_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// overlaps_C.h — C++ mapping of Fortran module "overlaps_C".
#pragma once
namespace overlaps_C {
constexpr int N_AO = 60;
// STO-6G coefficients and exponents, dimension(60,6), 1-based.
extern double ccc[N_AO + 1][7];
extern double zzz[N_AO + 1][7];
// STO-6G master tables (setupg.F90): allz/allc(6,6,2), 1-based [nqn][i][l].
extern double allz[7][7][3];
extern double allc[7][7][3];
// fact(0..17): factorials n!
extern double fact[18];
extern double cutof1, cutof2;
}
