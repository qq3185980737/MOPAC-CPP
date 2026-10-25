// ============================================================================
// 【中文导读】parameters_for_PM6_Sparkles_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// parameters_for_PM6_Sparkles_C.h — C++ translation of MOPAC 2016 "parameters_for_PM6_Sparkles_C.F90".
#pragma once

#include <vector>

namespace parameters_for_PM6_Sparkles_C {
extern std::vector<double> alp6sp;  // (107)
extern std::vector<double> gss6sp;  // (107)
extern std::vector<std::vector<double>> gues6sp1;  // (107,2)
extern std::vector<std::vector<double>> gues6sp2;  // (107,2)
extern std::vector<std::vector<double>> gues6sp3;  // (107,2)
}  // namespace parameters_for_PM6_Sparkles_C
