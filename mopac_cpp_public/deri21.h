// ============================================================================
// 【中文导读】deri21.h —— 解析导数/积分接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// deri21.h — C++ translation of MOPAC 2016 "deri21.F90".
#pragma once
// Least-squares (PCA) orthonormal basis extraction from set {A}.
// a: minear x nvar_nvo, column-major 0-based (F90 a(minear,nvar_nvo));
// vnert: nvar_nvo; pnert: nvar_nvo (sqrt of |eigenvalues|, ascending);
// b: minear x ncut, column-major 0-based output.
void deri21(const double* a, int nvar_nvo, int minear, int ifirst,
            double* vnert, double* pnert, double* b, int& ncut);
