// ============================================================================
// 【中文导读】ef.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// ef.h — C++ translation of MOPAC 2016 "ef.F90".
#pragma once
#include <vector>

// P-RFO/QA optimizer driver (ef.F90 "ef").
void ef(std::vector<double>& xparam, double& funct);
// Hessian update: iupd=1 Powell, iupd=2 BFGS (ef.F90 "updhes").
void updhes(std::vector<double>& svec, std::vector<double>& tvec,
            const std::vector<double>& grad, int nvar, int iupd,
            std::vector<double>& hess, const std::vector<double>& oldf,
            const std::vector<double>& d);
