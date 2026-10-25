// ============================================================================
// 【中文导读】prtpka.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// prtpka.h — C++ translation of MOPAC 2016 "prtpka.F90".
#pragma once

// Calculates and prints the pKa values for an organic molecule with an
// ionizable hydrogen attached to an oxygen atom.
void prtpka(int* ipKa_sorted, double* pKa_sorted, int* ipKa_unsorted,
            double* pKa_unsorted, int& no);

// pKa parameter set + NDDO parameter-table overrides (H/C/N/O/F/Cl/Br/I).
void Parameters_for_PKA(double& c1, double& c2, double& c3);
