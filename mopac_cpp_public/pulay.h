// ============================================================================
// 【中文导读】pulay.h —— SCF/Fock 接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// pulay.h — C++ translation.
#pragma once
void pulay(double* f, double* p, int n, double* fppf, double* fock,
           double* emat, int& lfock, int& nfock, int msize, bool& start,
           double& pl);
