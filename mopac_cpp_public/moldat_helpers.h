// ============================================================================
// 【中文导读】moldat_helpers.h —— MOZYME 接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// moldat_helpers.h — C++ translation of moldat.F90 internal subprograms.
#pragma once

// setcup: determine CUTOFP and unit-cell counts (moldat.F90 lines 897-1035).
void setcup();

// setup_nhco: MM correction to -(C=O)-(NH)- linkages (moldat.F90 1250-1323).
void setup_nhco(int& ii);
