// ============================================================================
// 【中文导读】parameters_for_mndod_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// parameters_for_mndod_C.h — C++ translation of MOPAC 2016 "parameters_for_mndod_C.F90".
#pragma once

#include <vector>

namespace parameters_for_mndod_C {
extern std::vector<double> alpd;  // (107)
extern std::vector<double> betadd;  // (107)
extern std::vector<double> betapd;  // (107)
extern std::vector<double> betasd;  // (107)
extern std::vector<double> gp2d;  // (107)
extern std::vector<double> gppd;  // (107)
extern std::vector<double> gspd;  // (107)
extern std::vector<double> gssd;  // (107)
extern std::vector<double> hspd;  // (107)
extern std::vector<double> poc_d;  // (107)
extern std::vector<double> uddd;  // (107)
extern std::vector<double> uppd;  // (107)
extern std::vector<double> ussd;  // (107)
extern std::vector<double> zdd;  // (107)
extern std::vector<double> zdnd;  // (107)
extern std::vector<double> zpd;  // (107)
extern std::vector<double> zpnd;  // (107)
extern std::vector<double> zsd;  // (107)
extern std::vector<double> zsnd;  // (107)
}  // namespace parameters_for_mndod_C
