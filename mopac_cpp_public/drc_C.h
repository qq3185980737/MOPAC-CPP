// ============================================================================
// 【中文导读】drc_C.h —— 反应路径/几何导数接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// drc_C.h — C++ mapping of Fortran module "drc_C" (subset used).
#pragma once
#include <vector>
namespace drc_C {
extern std::vector<double> vref;     // reference velocities (DRC)
extern std::vector<double> vref0;    // initial reference velocities
extern std::vector<double> allxyz;   // all cartesian coords along DRC
extern std::vector<double> allvel;   // all velocities along DRC
extern std::vector<double> xyz3;     // current xyz (3 x natoms)
extern std::vector<double> vel3;     // current velocities
extern std::vector<double> allgeo;   // all internal coords along DRC
extern std::vector<double> geo3;     // current internal coords
extern std::vector<double> parref;   // reference parameters
extern double time;          // DRC elapsed time (femtoseconds)
}
