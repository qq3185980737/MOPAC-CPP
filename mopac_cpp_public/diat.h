// ============================================================================
// 【中文导读】diat.h —— 解析导数/积分接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// diat.h — C++ translation of MOPAC 2016 "diat.F90".
#pragma once
#include <vector>

// Generalized Slater-type orbital overlap (pure); bfn external.
double ss_overlap(int na, int nb, int la, int lb, int m,
                  double ua, double ub, double r1, double a0);

// Diatomic overlap driver (skeleton).
void diat(int ni, int nj, const std::vector<double>& xj,
          std::vector<std::vector<double>>& di);
