// ============================================================================
// 【中文导读】analyt_C.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 analyt_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// analyt_C.cpp — storage for analyt_C.
#include "analyt_C.h"
namespace analyt_C {
int nztype[NZTYPE_MAX + 1] = {};
double ds[17] = {};
double dg[23] = {};
double dr[101] = {};
double g[23] = {};
double tx[4] = {}, ty[4] = {}, tz[4] = {};
double tdx[4] = {}, tdy[4] = {}, tdz[4] = {};
}
