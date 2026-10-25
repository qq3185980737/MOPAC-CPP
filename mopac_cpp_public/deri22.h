// ============================================================================
// 【中文导读】deri22.h —— 解析导数/积分接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// deri22.h — C++ translation of MOPAC 2016 "deri22.F90".
#pragma once
#include <vector>

// Build (D-A)*B supervector ab from density derivative b, plus the
// CI-active Fock diagonal blocks fci. b/ab/fci are column-major supervectors
// over the relaxation basis; bcol/abcol/fcicol select the 1-based column
// (F90 b(1,j) / ab(1,j) / fci(1,j)).
void deri22(const std::vector<std::vector<double>>& c,
            std::vector<double>& b, int bcol,
            std::vector<std::vector<double>>& work,
            std::vector<double>& foc2,
            std::vector<double>& ab, int abcol, int minear,
            std::vector<double>& fci, int fcicol,
            std::vector<double>& w,
            const std::vector<double>& diag,
            const std::vector<double>& scalar, int ninear);
