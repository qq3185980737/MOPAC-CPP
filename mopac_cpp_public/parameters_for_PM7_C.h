// ============================================================================
// 【中文导读】parameters_for_PM7_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// parameters_for_PM7_C.h — C++ translation of MOPAC 2016 "parameters_for_PM7_C.F90".
#pragma once

#include <vector>

namespace parameters_for_PM7_C {
extern std::vector<double> alp7;  // (107)
extern std::vector<double> betad7;  // (107)
extern std::vector<double> betap7;  // (107)
extern std::vector<double> betas7;  // (107)
extern std::vector<double> f0sd7;  // (107)
extern std::vector<double> g2sd7;  // (107)
extern std::vector<double> gp27;  // (107)
extern std::vector<double> gpp7;  // (107)
extern std::vector<double> gsp7;  // (107)
extern std::vector<double> gss7;  // (107)
extern std::vector<std::vector<double>> gues71;  // (107,4)
extern std::vector<std::vector<double>> gues72;  // (107,4)
extern std::vector<std::vector<double>> gues73;  // (107,4)
extern std::vector<double> hsp7;  // (107)
extern std::vector<double> poc_7;  // (107)
extern std::vector<double> polvo7;  // (107)
extern std::vector<double> udd7;  // (107)
extern std::vector<double> upp7;  // (107)
extern std::vector<double> uss7;  // (107)
extern std::vector<double> zd7;  // (107)
extern std::vector<double> zdn7;  // (107)
extern std::vector<double> zp7;  // (107)
extern std::vector<double> zpn7;  // (107)
extern std::vector<double> zs7;  // (107)
extern std::vector<double> zsn7;  // (107)
}  // namespace parameters_for_PM7_C
