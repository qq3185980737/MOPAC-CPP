// ============================================================================
// 【中文导读】mlmo.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// mlmo.h — C++ translation.
#pragma once
void mlmo(int& locc, int& lvir, int ii, int jj, int& nf_loc, int& ne,
          int& nocc, int& nvir, int* iz, int* ib, int* nce, int* ncf,
          int* ncocc, int* ncvir, const int* iorbs, int* icocc, int* icvir,
          double* cocc, double* cvir);
