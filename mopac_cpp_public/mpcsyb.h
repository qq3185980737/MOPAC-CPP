// ============================================================================
// 【中文导读】mpcsyb.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// mpcsyb.h — C++ translation.
#pragma once
#include <cstdio>
void mpcsyb(double* chr, int kchrge, double eionis, double& dip);

// F90 mpcpop(icok): Mulliken populations; isyb = open SYBYL output stream.
void mpcpop(int icok, FILE* isyb);
