// ============================================================================
// 【中文导读】eimp.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 eimp.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// eimp.cpp — C++ translation of MOPAC 2016 "eimp.F90".

#include "eimp.h"

#include "ijbo.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;


void eimp() {
    for (int i = 1; i <= numat; ++i)
        for (int j = 1; j < i; ++j) {
            int k = ijbo(i, j);
            if (k >= 0) {
                int l = MOZYME_C::iorbs[i] * MOZYME_C::iorbs[j];
                if (l != 0) {
                    l = k + l;
                    double sum = 0.0;
                    for (int m = k + 1; m <= l; ++m) sum += f[m] * f[m];
                    p[k + 1] = sum;
                }
            }
        }
}
