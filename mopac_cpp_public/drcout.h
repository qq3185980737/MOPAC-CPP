// ============================================================================
// 【中文导读】drcout.h —— 反应路径/几何导数接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// drcout.h — C++ translation of MOPAC 2016 "drcout.F90".
#pragma once
#include <string>
#include <vector>

void drcout(const std::vector<std::vector<double>>& xyz3,
            const std::vector<std::vector<double>>& geo3,
            const std::vector<std::vector<double>>& vel3,
            int nvar, double time,
            const std::vector<double>& escf3, const std::vector<double>& ekin3,
            const std::vector<double>& etot3, const std::vector<double>& xtot3,
            int iloop, const std::vector<double>& charge, double fract,
            const std::string& text1, const std::string& text2,
            int ii, int& jloop);
