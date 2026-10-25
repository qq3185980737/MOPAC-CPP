// ============================================================================
// 【中文导读】parameters_for_mndo_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// parameters_for_mndo_C.h — C++ translation of MOPAC 2016 "parameters_for_mndo_C.F90".
#pragma once

#include <vector>

namespace parameters_for_mndo_C {
extern std::vector<double> alpm;  // (107)
extern std::vector<double> betadm;  // (107)
extern std::vector<double> betapm;  // (107)
extern std::vector<double> betasm;  // (107)
extern std::vector<double> f0sdm;  // (107)
extern std::vector<double> g2sdm;  // (107)
extern std::vector<double> gp2m;  // (107)
extern std::vector<double> gppm;  // (107)
extern std::vector<double> gspm;  // (107)
extern std::vector<double> gssm;  // (107)
extern std::vector<std::vector<double>> guesm1;  // (107,4)
extern std::vector<std::vector<double>> guesm2;  // (107,4)
extern std::vector<std::vector<double>> guesm3;  // (107,4)
extern std::vector<double> hspm;  // (107)
extern std::vector<double> pocm;  // (107)
extern std::vector<double> polvolm;  // (107)
extern std::vector<double> polvom;  // (107)
extern std::vector<double> uddm;  // (107)
extern std::vector<double> uppm;  // (107)
extern std::vector<double> ussm;  // (107)
extern std::vector<double> zdm;  // (107)
extern std::vector<double> zdnm;  // (107)
extern std::vector<double> zpm;  // (107)
extern std::vector<double> zpnm;  // (107)
extern std::vector<double> zsm;  // (107)
extern std::vector<double> zsnm;  // (107)
}  // namespace parameters_for_mndo_C
