// ============================================================================
// 【中文导读】ef_C.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 ef_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

#include "ef_C.h"
namespace ef_C {
int nstep=0, negreq=0, iprnt=0, ef_mode=0;
double ddx=1.0, xlamd=0, xlamd0=0, skal=1, rmin=0, rmax=0, omin=0;
double x0=0, x1=0, x2=0; int iloop=0;
std::vector<std::vector<double>> alparm;
std::vector<std::vector<double>> hess;
std::vector<std::vector<double>> bmat;
std::vector<double> u, oldhss, oldu;
std::vector<std::vector<double>> pmat, uc, oldf, d, vmode;
std::vector<double> hessc;
}
