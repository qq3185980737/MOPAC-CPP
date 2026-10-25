// ============================================================================
// 【中文导读】helecz.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 helecz.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// helecz.cpp
#include "helecz.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
using namespace molkst_C;
using namespace MOZYME_C;
using namespace common_arrays_C;
extern int ijbo(int, int);
extern "C" int ijbo_(int* i, int* j) { return ijbo(*i, *j); }
double helecz() {
    double ed = 0.0, ee = 0.0;
    for (int i = 1; i <= numat; ++i) {
        for (int j = 1; j <= i - 1; ++j) {
            int k = ijbo_(&i, &j);
            if (k >= 0) {
                int l = k + iorbs[i] * iorbs[j];
                for (int m = k + 1; m <= l; ++m)
                    ee += p[m] * (h[m] + f[m]);
            }
        }
        int k = ijbo_(&i, &i);
        for (int l1 = 1; l1 <= iorbs[i]; ++l1) {
            for (int l2 = 1; l2 <= l1 - 1; ++l2) {
                k++;
                ee += p[k] * (h[k] + f[k]);
            }
            k++;
            ed += p[k] * (h[k] + f[k]);
        }
    }
    ee += 0.5 * ed;
    return ee;
}
