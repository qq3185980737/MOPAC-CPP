// ============================================================================
// 【中文导读】aabacd.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// aabacd.h — C++ translation of MOPAC 2016 "aabacd.F90".
// Two microstates differing by two alpha MOs; xy flat column-major 4-D,
// F90 1-based semantics.
#pragma once
double aabacd(const int* iocca1, const int* ioccb1, const int* iocca2,
              const int* ioccb2, int nmos, const double* xy);
