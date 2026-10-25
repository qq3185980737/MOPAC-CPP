// ============================================================================
// 【中文导读】mecih.cpp —— CI/激发态与极化
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mecih.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mecih.cpp — build the packed canonical CI matrix CIMAT.
#include "mecih.h"

#include <cmath>
#include <vector>

#include "aabbcd.h"
#include "aababc.h"
#include "aabacd.h"
#include "babbcd.h"
#include "babbbc.h"
#include "meci_C.h"

using namespace meci_C;

void mecih(double* diag, double* cimat, int nmos, const int* /*lab_in*/,
           double* xy) {
    int ik = 0;
    for (iiloop = 1; iiloop <= lab; ++iiloop) {
        is = 2;
        for (jloop = 1; jloop <= iiloop; ++jloop) {
            ++ik;
            cimat[ik - 1] = 0.0;   // physical 0-based packed slot (rsp reads a[0]=slot1)
            int ix = 0, iy = 0;
            for (int j = 1; j <= nmos; ++j) {
                ix += std::abs(microa[j][iiloop] - microa[j][jloop]);
                iy += std::abs(microb[j][iiloop] - microb[j][jloop]);
            }
            if (ix + iy > 4 || nalmat[iiloop] != nalmat[jloop]) {
                continue;
            }

            std::vector<int> ma1(nmos + 1), mb1(nmos + 1), ma2(nmos + 1),
                mb2(nmos + 1);
            for (int j = 1; j <= nmos; ++j) {
                ma1[j] = microa[j][iiloop];
                mb1[j] = microb[j][iiloop];
                ma2[j] = microa[j][jloop];
                mb2[j] = microb[j][jloop];
            }

            if (ix + iy == 4) {
                if (ix == 0)
                    cimat[ik - 1] = babbcd(&ma1[0], &mb1[0], &ma2[0], &mb2[0],
                                           nmos, xy);
                else if (ix == 2)
                    cimat[ik - 1] = aabbcd(&ma1[0], &mb1[0], &ma2[0], &mb2[0],
                                           nmos, xy);
                else
                    cimat[ik - 1] = aabacd(&ma1[0], &mb1[0], &ma2[0], &mb2[0],
                                           nmos, xy);
            } else if (ix == 2) {
                cimat[ik - 1] = aababc(&ma1[0], &mb1[0], &ma2[0], nmos, xy);
            } else if (iy == 2) {
                cimat[ik - 1] = babbbc(&ma1[0], &mb1[0], &mb2[0], nmos, xy);
            } else {
                cimat[ik - 1] = diag[iiloop];
            }
        }
        ispqr[iiloop][1] = is - 1;
    }
}
