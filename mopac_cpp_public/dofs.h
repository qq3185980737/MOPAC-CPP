// ============================================================================
// 【中文导读】dofs.h —— 解析导数/积分接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// dofs.h — C++ translation of MOPAC 2016 "dofs.F90".
#pragma once
#include <vector>

// Density-of-states histogram into dd(m). eref(mono3,n).
void dofs(std::vector<std::vector<double>>& eref, int mono3, int n,
          std::vector<double>& dd, int m, double bottom, double top);
