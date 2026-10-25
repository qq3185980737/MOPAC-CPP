// ============================================================================
// 【中文导读】prttim.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 prttim.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// prttim.cpp
#include "prttim.h"
void prttim(double tleft, double& tprt, char& txt) {
    tprt=tleft; txt='S';
    if (tprt>=604800.0) { tprt/=604800.0; txt='W'; }
    else if (tprt>=86400.0) { tprt/=86400.0; txt='D'; }
    else if (tprt>=3600.0) { tprt/=3600.0; txt='H'; }
    else if (tprt>=60.0) { tprt/=60.0; txt='M'; }
}
