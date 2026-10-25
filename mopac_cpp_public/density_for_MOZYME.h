// ============================================================================
// 【中文导读】density_for_MOZYME.h —— 密度矩阵接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// density_for_MOZYME.h — C++ translation of MOPAC 2016.
#pragma once
#include <vector>

// Build packed density p from compressed occupied LMOs.
void density_for_MOZYME(std::vector<double>& p, int mode, int nclose_loc,
                        const std::vector<double>& partpin);
