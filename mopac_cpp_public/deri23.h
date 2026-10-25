// ============================================================================
// 【中文导读】deri23.h —— 解析导数/积分接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// deri23.h — C++ translation of MOPAC 2016 "deri23.F90".
#pragma once
#include <vector>

// Unpack CI-active MO derivatives into cmo; eigenvalue relaxation into emo.
// f/fd/fci are column-major supervectors over the full basis (1-based
// padding); ivar selects the current geometric-variable column so the caller
// passes the same container for every column (as F90 f(1,ivar)).
void deri23(const std::vector<double>& f,
            const std::vector<double>& fd,
            const std::vector<double>& e,
            const std::vector<double>& fci,
            std::vector<std::vector<double>>& cmo,
            std::vector<double>& emo, int minear, int ninear, int ivar);
