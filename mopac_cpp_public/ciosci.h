// ============================================================================
// 【中文导读】ciosci.h —— CI/激发态接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// ciosci.h — C++ translation of MOPAC 2016 "ciosci.F90".
#pragma once
#include <vector>

// Oscillator-strength / dipole matrix between CI states.
// vects(norbs,norbs), oscil(3,lab), conf(lab*lab).
void ciosci(const std::vector<std::vector<double>>& vects, int lroot,
            std::vector<std::vector<double>>& oscil,
            const std::vector<double>& conf_in);
