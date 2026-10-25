// ============================================================================
// 【中文导读】dfield.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 dfield.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// dfield.cpp — C++ translation of MOPAC 2016 "dfield.F90".
#include "dfield.h"

#include <vector>

#include "common_arrays_C.h"
#include "chrge.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace parameters_C;

void dfield() {
    std::vector<double> q2(numat + 1, 0.0);
    chrge(p, q2);
    for (int i = 1; i <= numat; ++i) q2[i] = tore[nat[i]] - q2[i];
    double fldcon = ev / a0 * fpc_9;
    for (int i = 1; i <= numat; ++i) {
        dxyz[3 * (i - 1) + 1] += efield[1] * q2[i] * fldcon;
        dxyz[3 * (i - 1) + 2] += efield[2] * q2[i] * fldcon;
        dxyz[3 * (i - 1) + 3] += efield[3] * q2[i] * fldcon;
    }
}
