// ============================================================================
// 【中文导读】extvdw_for_MOZYME.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// extvdw_for_MOZYME.h — C++ translation of "extvdw_for_MOZYME" (from
// geochk.F90, lines 1490-1579).  Fills per-atom van der Waals radii,
// applying keyword "METAL" and "VDWM(:sym=value:...)" overrides.
#pragma once
#include <vector>

void extvdw_for_MOZYME(std::vector<double>& radius,
                       const std::vector<double>& refvdw);
