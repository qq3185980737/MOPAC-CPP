// ============================================================================
// 【中文导读】rotate_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// rotate_C.h — C++ translation of MOPAC 2016 "rotate_C.F90".
// Fortran: ccore(4,2) with equivalence ccore(i,1) <-> css_i (i=1..4) and
// ccore(i,2) <-> css2_i. Column dimension kept 1-based -> [3].
#pragma once

namespace rotate_C {
// thread_local: rotate() runs concurrently under the dcart OpenMP atom-pair
// loop; these Fortran-equivalence scratch globals must be per-thread.
extern thread_local double css1, csp1, cpps1, cppp1;  // first  interaction (s, sp, pps, ppp)
extern thread_local double css2, csp2, cpps2, cppp2;  // second interaction
extern thread_local double ccore[4][3];
}  // namespace rotate_C
