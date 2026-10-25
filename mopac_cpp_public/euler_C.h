// ============================================================================
// 【中文导读】euler_C.h —— 数值/几何工具接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// euler_C.h — C++ translation of MOPAC 2016 "euler_C.F90".
//
// Fortran module euler_C holds: id, l1l..l3l, l1u..l3u, tvec(3,3).
// In this C++ port those globals are split by existing architecture:
//   id, l1u, l2u, l3u  -> molkst_C  (molkst_C.h)
//   tvec(3,3)          -> common_arrays_C (vector form, periodic lattice)
//   l1l, l2l, l3l      -> this namespace (the only members defined here)
#pragma once

namespace euler_C {
extern int l1l, l2l, l3l;
}  // namespace euler_C
