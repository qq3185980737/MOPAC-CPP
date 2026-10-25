// ============================================================================
// 【中文导读】freqcy.h —— 力场/振动/热力学接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// freqcy.h — C++ translation of MOPAC 2016 "freqcy.F90".
#pragma once
#include <vector>

// Mass-weight lower-triangle packed Hessian (1-based). fmatrx modified in place.
// oldf (same size) receives the unweighted copy.
void freqcy_mass_weight(std::vector<double>& fmatrx,
                        std::vector<double>& oldf,
                        const std::vector<double>& wtmass, int nvar);
void freqcy(std::vector<double>& fmatrx, std::vector<double>& freq,
            std::vector<double>& travel, bool eorc,
            std::vector<std::vector<double>>& deldip,
            std::vector<double>& ff, std::vector<double>& oldf, bool ts);
