// ============================================================================
// 【中文导读】axis.h —— 力场/振动/热力学接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// axis.h — C++ translation of MOPAC 2016 "axis.F90".
#pragma once

// axis: compute the three principal moments of inertia and molecular weight.
// Rotational constants (cm^-1) are returned in a, b, c; eigenvectors evec(1..3,1..3).
void axis(double& a, double& b, double& c, double evec[4][4]);
