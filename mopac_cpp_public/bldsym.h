// ============================================================================
// 【中文导读】bldsym.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// bldsym.h — C++ translation of MOPAC 2016 "bldsym.F90".
#pragma once

// bldsym: build the 3x3 point-group operation matrix for ioper into
// symmetry_C::elem(.., .., iplace).
void bldsym(int ioper, int iplace);
