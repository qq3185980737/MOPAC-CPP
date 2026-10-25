// ============================================================================
// 【中文导读】setupk.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 setupk.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// setupk.cpp — C++ translation of MOPAC 2016 "setupk.F90".
// Uses MOZYME_C (shared) arrays icocc/ncf/nncf/kopt (1-based).
#include "setupk.h"
#include <algorithm>
#include "molkst_C.h"
#include "MOZYME_C.h"
using molkst_C::numat;
using namespace MOZYME_C;
void setupk(int nocc1) {
    std::fill(kopt.begin(), kopt.end(), 0);
    for (int i = 1; i <= nocc1; ++i) {
        int j = nncf[i];
        for (int k = 1; k <= ncf[i]; ++k)
            kopt[icocc[k + j]] = 1;
    }
    int l = 0;
    for (int i = 1; i <= numat; ++i) {
        if (kopt[i] == 1) { ++l; kopt[l] = i; }
    }
    if (l != numat) kopt[l + 1] = 0;
}
