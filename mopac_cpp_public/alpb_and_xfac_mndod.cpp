// ============================================================================
// 【中文导读】alpb_and_xfac_mndod.cpp —— ALPB/xfac 参数表（按方法/元素）
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 alpb_and_xfac_mndod.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// alpb_and_xfac_mndod.cpp — C++ translation of MOPAC 2016 "alpb_and_xfac_mndod.F90" (alpb_and_xfac_mndod).
#include "alpb_and_xfac.h"
#include "parameters_C.h"

void alpb_and_xfac_mndod() {
    using parameters_C::alpb;
    using parameters_C::xfac;
    alpb[11][1] = 1.05225212e0;
    xfac[11][1] = 1.00000000e0;
    alpb[11][6] = 1.05225212e0;
    xfac[11][6] = 1.00000000e0;
    alpb[12][1] = 1.35052992e0;
    xfac[12][1] = 1.00000000e0;
    alpb[12][6] = 1.48172071e0;
    xfac[12][6] = 1.00000000e0;
    alpb[16][12] = 1.48172071e0;
    xfac[16][12] = 1.00000000e0;
    alpb[13][1] = 1.38788000e0;
    xfac[13][1] = 1.00000000e0;
    alpb[13][6] = 1.38788000e0;
    xfac[13][6] = 1.00000000e0;
    alpb[13][13] = 1.38788000e0;
    xfac[13][13] = 1.00000000e0;
}
