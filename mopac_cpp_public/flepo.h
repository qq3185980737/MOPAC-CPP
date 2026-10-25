// ============================================================================
// 【中文导读】flepo.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// flepo.h — C++ translation of MOPAC 2016 "flepo.F90".
#pragma once
#include <vector>

// Inverse-Hessian update (packed lower triangle, length nvar*(nvar+1)/2).
// dfp=true: DFP; false: BFGS. hesinv updated in place.
void flepo_update_H(int nvar, const std::vector<double>& xvar,
                    const std::vector<double>& gvar,
                    const std::vector<double>& gg,
                    double sy, double yhy, bool dfp,
                    std::vector<double>& hesinv);
void flepo(std::vector<double>& xparam, int nvar, double& funct1);
