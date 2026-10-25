// ============================================================================
// 【中文导读】hbonds.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 hbonds.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// hbonds.cpp
#include "hbonds.h"
#include <algorithm>
#include <vector>
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "ijbo.h"

using namespace molkst_C;
using namespace MOZYME_C;
using namespace common_arrays_C;

void hbonds(const double* fao, int nocc, int nvir, int* iused, int& nij_loc, double cutoff) {
    nij_loc = 0;
    for (int i = 1; i <= numat; ++i) {
        for (int j = 1; j <= i - 1; ++j) {
            if (ijbo(i, j) >= 0) {
                int ll = ijbo(i, j);
                double sumf = 0.0, sump = 0.0;
                for (int k = 1; k <= iorbs[i]; ++k) {
                    for (int l = 1; l <= iorbs[j]; ++l) {
                        ++ll;
                        sumf += fao[ll] * fao[ll];
                        sump += p[ll] * p[ll];
                    }
                }
                if (sump < 1e-10 && sumf > cutoff) {
                    int m = 0;
                    for (int k = 1; k <= nocc; ++k) {
                        int kk = std::min(ncf[k] - 1, 1);
                        int l = nncf[k] + 1;
                        if (icocc[l] == i || icocc[l + kk] == i ||
                            icocc[l] == j || icocc[l + kk] == j) {
                            ++m;
                            iused[m] = k;
                        }
                    }
                    for (int k = 1; k <= nvir; ++k) {
                        int kk = std::min(nce[k] - 1, 1);
                        int l = nnce[k] + 1;
                        if (icvir[l] == i || icvir[l + kk] == i ||
                            icvir[l] == j || icvir[l + kk] == j) {
                            for (int ll2 = 1; ll2 <= m; ++ll2) {
                                ++nij_loc;
                                fmo[nij_loc] = 0.1;
                                ifmo[2][nij_loc] = iused[ll2];
                                ifmo[1][nij_loc] = k;
                            }
                        }
                    }
                }
            }
        }
    }
}
