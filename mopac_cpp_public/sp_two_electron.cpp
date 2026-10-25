// ============================================================================
// 【中文导读】sp_two_electron.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 sp_two_electron.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// sp_two_electron.cpp
#include "sp_two_electron.h"
#include "parameters_C.h"
#include "mndod_C.h"
using namespace parameters_C;
using mndod_C::iii;
extern double rsc(int, int, double, int, double, int, double, int, double);
void sp_two_electron() {
    for (int ni=1; ni<=80; ++ni) {
        int ns = iii[ni];
        double es = zsn[ni], ep = zpn[ni];
        if (es<1e-4 || ep<1e-4 || main_group[ni]) continue;
        gss[ni] = rsc(0, ns, es, ns, es, ns, es, ns, es);
        gsp[ni] = rsc(0, ns, es, ns, es, ns, ep, ns, ep);
        hsp[ni] = rsc(1, ns, es, ns, ep, ns, es, ns, ep)/3.0;
        double r033 = rsc(0, ns, ep, ns, ep, ns, ep, ns, ep);
        double r233 = rsc(2, ns, ep, ns, ep, ns, ep, ns, ep);
        gpp[ni] = r033 + 0.16*r233;
        gp2[ni] = r033 - 0.08*r233;
    }
}
