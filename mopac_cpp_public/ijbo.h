// ============================================================================
// 【中文导读】ijbo.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// ijbo.h — C++ translation of ijbo.F90 (MOPAC 2016).
#pragma once

// Given atom numbers ii and jj:
//   IJBO = -1 if the atoms are separated by more than SQRT(CUTOF1), else
//   IJBO = -2 if the atoms are separated by more than SQRT(CUTOF2), else
//   IJBO = N, where "N+1" is the address of the starting location in the
//           H, P, and F arrays.
int ijbo(int ii, int jj);