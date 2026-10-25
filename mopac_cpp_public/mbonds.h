// ============================================================================
// 【中文导读】mbonds.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// mbonds.h — C++ translation.
#pragma once
void mbonds(int locc, int lvir, const double* f, double* catom,
            const int* nfirst, int ii, int jj, bool* ui, bool* uj,
            bool& lok, const int* iorbs, double* cocc, double* cvir,
            int cocc_dim, int cvir_dim, int numat, int norbs, int morb, int mpack1);
