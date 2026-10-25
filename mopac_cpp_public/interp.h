// ============================================================================
// 【中文导读】interp.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// interp.h — C++ translation.
#pragma once
void interp(int np, int nq, int& mode, double e, double* fp, double* cp,
            double* theta, double* vec_interp, double* fock_interp,
            double* p_interp, double* h_interp, double* vecl, double& eold_ref);
void spline(double* x, double* f, double* df, double xhigh, double xlow,
            double& xmin, int npnts);
