// ============================================================================
// 【中文导读】switch.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// switch.h — C++ translation of MOPAC 2016 "switch.F90".
#pragma once

// Copies reference parameter sets (parameters_for_*_C) into the runtime
// arrays of parameters_C according to the active method flags, then runs
// the alpb_and_xfac_<model> fill, symmetrization and fractional-metal fixes.
void switch_method();
