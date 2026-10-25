// ============================================================================
// 【中文导读】compfg.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// compfg.h — C++ translation of MOPAC 2016 "compfg.F90".
// The SCF driver is a skeleton; xfac_value (core-core energy) is fully implemented.
#pragma once
#include <vector>

void compfg(const std::vector<double>& xparam, bool int_flag, double& escf,
            bool fulscf, std::vector<double>& grad, bool lgrad);

// Core-core repulsion used by xfac keyword.
double xfac_value();
