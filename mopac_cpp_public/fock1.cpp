// ============================================================================
// 【中文导读】fock1.cpp —— 单中心 Fock 矩阵项（fock1.F90）
// ----------------------------------------------------------------------------
// 由密度矩阵 p 与单中心双电子积分（w/wk 表）累加 Fock 矩阵 f 的单原子贡献，
// 是 buildf 的补充路径之一。1-based 打包下三角：f/p 索引 i(i-1)/2+j，
// 参数表来自 parameters_C。
// ============================================================================

// fock1.cpp — C++ translation of MOPAC 2016 "fock1.F90".

#include "fock1.h"

#include "common_arrays_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace parameters_C;

void fock1(std::vector<double>& Fmat, const std::vector<double>& ptot,
           const std::vector<double>& PAa, const std::vector<double>& PBb) {
    for (int ii = 1; ii <= numat; ++ii) {
        int ia = nfirst[ii], ib = nlast[ii], ni = nat[ii];
        double ptpop = 0.0, papop = 0.0;
        int nbases = ib - ia + 1;
        if (nbases >= 3) {
            ptpop = ptot[(ib * (ib + 1)) / 2] + ptot[((ib - 1) * ib) / 2] +
                    ptot[((ib - 2) * (ib - 1)) / 2];
            papop = PAa[(ib * (ib + 1)) / 2] + PAa[((ib - 1) * ib) / 2] +
                    PAa[((ib - 2) * (ib - 1)) / 2];
        }
        int ka = (ia * (ia + 1)) / 2;
        Fmat[ka] += PBb[ka] * gss[ni] + ptpop * gsp[ni] - papop * hsp[ni];
        if (nbases != 1) {
            int l = ka;
            for (int j = ia + 1; j <= ib; ++j) {
                int m = l + ia;
                l = l + j;
                Fmat[l] += ptot[ka] * gsp[ni] - PAa[ka] * hsp[ni] + PBb[l] * gpp[ni] +
                        (ptpop - ptot[l]) * gp2[ni] -
                        0.5 * (papop - PAa[l]) * (gpp[ni] - gp2[ni]);
                Fmat[m] += 2.0 * ptot[m] * hsp[ni] - PAa[m] * (hsp[ni] + gsp[ni]);
            }
            for (int j = ia + 1; j <= ib - 1; ++j)
                for (int l2 = j + 1; l2 <= ib; ++l2) {
                    int m = (l2 * (l2 - 1)) / 2 + j;
                    Fmat[m] += ptot[m] * (gpp[ni] - gp2[ni]) -
                            0.5 * PAa[m] * (gpp[ni] + gp2[ni]);
                }
        }
    }
}
