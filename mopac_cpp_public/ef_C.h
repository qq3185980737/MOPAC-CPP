// ============================================================================
// 【中文导读】ef_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// ef_C.h — globals from ef.C (extended for setup_mopac_arrays).
#pragma once
#include <vector>
namespace ef_C {
extern int nstep, negreq, iprnt, ef_mode;
extern double ddx, xlamd, xlamd0, skal, rmin, rmax, omin; extern double x0, x1, x2; extern int iloop; extern std::vector<std::vector<double>> alparm;
// EF optimizer workspaces (setup_mopac_arrays).
extern std::vector<std::vector<double>> hess;   // Hessian (3*nvar x 3*nvar)
extern std::vector<std::vector<double>> bmat;   // ef B-matrix (alias of ef_C bmat)
extern std::vector<double> u;                   // unitary transform scratch
extern std::vector<double> oldhss, oldu;        // old Hessian / old u
extern std::vector<std::vector<double>> pmat;   // p matrix
extern std::vector<std::vector<double>> uc;     // u * c
extern std::vector<double> hessc;               // Hessian * c
extern std::vector<std::vector<double>> oldf;   // old f matrix
extern std::vector<std::vector<double>> d;      // d matrix
extern std::vector<std::vector<double>> vmode;  // vibrational modes
}
