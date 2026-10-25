// ============================================================================
// 【中文导读】parameters_for_PM6_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// parameters_for_PM6_C.h — C++ translation of MOPAC 2016 "parameters_for_PM6_C.F90".
#pragma once

#include <vector>

namespace parameters_for_PM6_C {
extern std::vector<double> alp6;  // (107)
extern std::vector<double> betad6;  // (107)
extern std::vector<double> betap6;  // (107)
extern std::vector<double> betas6;  // (107)
extern std::vector<double> f0sd6;  // (107)
extern std::vector<double> g2sd6;  // (107)
extern std::vector<double> gp26;  // (107)
extern std::vector<double> gpp6;  // (107)
extern std::vector<double> gsp6;  // (107)
extern std::vector<double> gss6;  // (107)
extern std::vector<std::vector<double>> gues61;  // (107,4)
extern std::vector<std::vector<double>> gues62;  // (107,4)
extern std::vector<std::vector<double>> gues63;  // (107,4)
extern std::vector<double> hsp6;  // (107)
extern std::vector<double> poc_6;  // (107)
extern std::vector<double> polvo6;  // (107)
extern std::vector<double> udd6;  // (107)
extern std::vector<double> upp6;  // (107)
extern std::vector<double> uss6;  // (107)
extern std::vector<double> zd6;  // (107)
extern std::vector<double> zdn6;  // (107)
extern std::vector<double> zp6;  // (107)
extern std::vector<double> zpn6;  // (107)
extern std::vector<double> zs6;  // (107)
extern std::vector<double> zsn6;  // (107)
}  // namespace parameters_for_PM6_C
