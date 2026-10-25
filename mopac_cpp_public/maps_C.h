// ============================================================================
// 【中文导读】maps_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

#pragma once
#include <vector>
namespace maps_C {
extern double rxn_coord, rc_escf, ekin; extern int lparam, latom; extern std::vector<double> react; extern int kloop;
extern int lpara1, latom1, lpara2, latom2;   // reaction-path marker atoms (geout.F90)
// grid.F90 state.
extern double rxn_coord1, rxn_coord2;
extern int ione, ijlp, ilp, jlp, jlp1;
extern std::vector<double> surf;   // grid energies per point
}
