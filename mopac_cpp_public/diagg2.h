// ============================================================================
// 【中文导读】diagg2.h —— 对角化/本征值接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// diagg2.h — C++ translation of MOPAC 2016 "diagg2.F90".
#pragma once

// Jacobi annihilation of occupied-virtual LMO couplings with adaptive LMO
// expansion.  Pointer-based signature matches addhb.cpp's forward
// declaration: eigv is 0-based over nvirtual eigenvalues; storei/storej are
// (norbs) scratch; iused/latoms are (numat).
void diagg2(int nocc, int nvir, const double* eigv, int* iused,
            char* latoms, int nij, int idiagg, double* storei,
            double* storej);
