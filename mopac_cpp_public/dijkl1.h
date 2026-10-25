// ============================================================================
// 【中文导读】dijkl1.h —— 解析导数/积分接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// dijkl1.h — C++ translation of MOPAC 2016 "dijkl1.F90".
#pragma once
#include <vector>

// Two-electron integral derivatives over active MOs w.r.t. atom nati.
void dijkl1(const std::vector<std::vector<double>>& c, int n, int nati,
            const std::vector<double>& w,
            std::vector<double>& cij, std::vector<double>& wcij,
            std::vector<double>& ckl,
            std::vector<double>& xy);
