// ============================================================================
// 【中文导读】force.h —— 力场/振动/热力学接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// force.h — C++ translation of MOPAC 2016 "force.F90".
#pragma once
#include <vector>

// phase_lock: force the phase of each eigenvector so the largest
// coefficient is positive. vecs is packed n*n column-major (1-based).
void phase_lock(std::vector<double>& vecs, int n);
void force();
