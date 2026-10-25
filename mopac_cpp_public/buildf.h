// ============================================================================
// 【中文导读】buildf.h —— Fock/CI 超向量接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// buildf.h — C++ translation of MOPAC 2016 "buildf.F90".
#pragma once
#include <vector>

// buildf: combine a partial Fock matrix with the one-electron matrix h,
// then hand off to fock2z.
//   mode=-1: f = partf - h
//   mode= 0: f = h
//   mode= 1: f = partf + h
void buildf(std::vector<double>& f, const std::vector<double>& partf, int mode);
