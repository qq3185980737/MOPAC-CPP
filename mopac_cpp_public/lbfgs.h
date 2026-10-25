// ============================================================================
// 【中文导读】lbfgs.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// lbfgs.h — C++ translation of MOPAC 2016 "lbfgs.F90".
#pragma once
#include <string>
#include <vector>

// Driver used by run_mopac dispatch ("LBFGS"). xparam is 1-based (size nvar+1).
void lbfgs(double* xparam, double& escf);

// L-BFGS-B restart-file save (mode=1) / restore (mode=0).
void lbfsav(double tt0, int mode, std::vector<double>& wa, int nwa,
            std::vector<int>& iwa, int niwa, std::string& task,
            std::string& csave, std::vector<int>& lsave, std::vector<int>& isave,
            std::vector<double>& dsave, int& nstep, double& escf);

// L-BFGS-B entry point (called by lbfgs_TS). All arrays are 1-based.
void setulb(int n, int m, double* x, double* l, double* u, int* nbd,
            double& f, double* g, double factr, double pgtol,
            std::vector<double>& wa, std::vector<int>& iwa, std::string& task,
            int iprint, std::string& csave, std::vector<int>& lsave,
            std::vector<int>& isave, std::vector<double>& dsave);
