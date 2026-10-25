// ============================================================================
// 【中文导读】parameters_for_PM3_Sparkles_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// parameters_for_PM3_Sparkles_C.h — C++ translation of MOPAC 2016 "parameters_for_PM3_Sparkles_C.F90".
#pragma once

#include <vector>

namespace parameters_for_PM3_Sparkles_C {
extern std::vector<double> alpPM3sp;  // (107)
extern std::vector<double> gssPM3sp;  // (107)
extern std::vector<std::vector<double>> guesPM3sp1;  // (107,4)
extern std::vector<std::vector<double>> guesPM3sp2;  // (107,4)
extern std::vector<std::vector<double>> guesPM3sp3;  // (107,4)
}  // namespace parameters_for_PM3_Sparkles_C
