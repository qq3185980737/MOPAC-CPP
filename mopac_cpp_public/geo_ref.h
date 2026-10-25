// ============================================================================
// 【中文导读】geo_ref.h —— 输入/几何接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// geo_ref.h — C++ translation of geo_ref.F90 (MOPAC 2016).
#pragma once

// geo_ref reads in the geometry reference file and makes it as similar as
// possible to the current geometry: PDB labels are re-sequenced, the
// geometries are rotated/translated into maximum coincidence (RMS minimum),
// and pairs of atoms are swapped when that lowers the RMS difference.
void geo_ref();
