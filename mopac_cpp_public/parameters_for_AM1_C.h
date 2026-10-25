// ============================================================================
// 【中文导读】parameters_for_AM1_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// parameters_for_AM1_C.h — C++ translation of MOPAC 2016 "parameters_for_AM1_C.F90".
#pragma once

#include <vector>

namespace parameters_for_AM1_C {
extern std::vector<double> alpam1;  // (107)
extern std::vector<double> betada;  // (107)
extern std::vector<double> betapa;  // (107)
extern std::vector<double> betasa;  // (107)
extern std::vector<double> f0sdam1;  // (107)
extern std::vector<double> g2sdam1;  // (107)
extern std::vector<double> gp2am1;  // (107)
extern std::vector<double> gppam1;  // (107)
extern std::vector<double> gspam1;  // (107)
extern std::vector<double> gssam1;  // (107)
extern std::vector<std::vector<double>> guesa1;  // (107,4)
extern std::vector<std::vector<double>> guesa2;  // (107,4)
extern std::vector<std::vector<double>> guesa3;  // (107,4)
extern std::vector<double> hspam1;  // (107)
extern std::vector<double> polvolam1;  // (107)
extern std::vector<double> uddam1;  // (107)
extern std::vector<double> uppam1;  // (107)
extern std::vector<double> ussam1;  // (107)
extern std::vector<double> zdam1;  // (107)
extern std::vector<double> zdnam1;  // (107)
extern std::vector<double> zpam1;  // (107)
extern std::vector<double> zpnam1;  // (107)
extern std::vector<double> zsam1;  // (107)
extern std::vector<double> zsnam1;  // (107)
}  // namespace parameters_for_AM1_C
