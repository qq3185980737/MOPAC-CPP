// ============================================================================
// 【中文导读】deri1.h —— 解析导数/积分接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// deri1.h — C++ translation of MOPAC 2016 "deri1.F90".
#pragma once
#include <vector>

// Non-relaxed derivative wrt one Cartesian coordinate (skeleton).
void deri1(int number, double& grad,
           std::vector<double>& fv, int minear,
           std::vector<double>& fd,
           const std::vector<double>& scalar,
           std::vector<std::vector<double>>& work);
