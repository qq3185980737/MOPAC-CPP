// ============================================================================
// 【中文导读】lyse.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 lyse.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// lyse.cpp
#include "lyse.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
using namespace common_arrays_C;
using namespace molkst_C;
extern double distance(int, int);
extern "C" double distance_(int* a, int* b) { return distance(*a, *b); }
void lyse() {
    for (int i = 1; i <= numat; ++i) {
        if (nat[i] != 1) {
            int l = 0;
            if (nat[i] == 6) {
                for (int k = 1; k <= nbonds[i]; ++k) { l++; ibonds[l][i] = ibonds[k][i]; }
            } else {
                bool special = false;
                if (nat[i] == 8) {
                    for (int k = 1; k <= nbonds[i]; ++k) {
                        int nbr = ibonds[k][i];
                        if ((nat[nbr] == 16 || nat[nbr] == 15) && nbonds[nbr] == 4) special = true;
                    }
                }
                if (special) { l = nbonds[i]; goto done_i; }
                l = 0;
                for (int k = 1; k <= nbonds[i]; ++k) {
                    if (nat[ibonds[k][i]] != 16) { l++; ibonds[l][i] = ibonds[k][i]; }
                }
                {
                    int ll = 0;
                    for (int k = 1; k <= l; ++k) {
                        int j = ibonds[k][i];
                        if (nat[j] == 7 || nat[j] == 6) ll = 1;
                    }
                    if (ll == 0) continue;  // F90: cycle skips nbonds(i)=l
                    int ll2 = l; l = 0;
                    for (int k = 1; k <= ll2; ++k) {
                        int j = ibonds[k][i];
                        int n_j = nat[j];
                        if (n_j == 1 || n_j == 6 || n_j == 7 || n_j == 8 || n_j == 16) {
                            l++; ibonds[l][i] = ibonds[k][i];
                        } else {
                            int m = 0;
                            for (int n = 1; n <= nbonds[j]; ++n) {
                                if (ibonds[n][j] != i) { m++; ibonds[m][j] = ibonds[n][j]; }
                            }
                            nbonds[j]--;
                        }
                    }
                }
            }
done_i:
            nbonds[i] = l;
        } else {
            if (nbonds[i] > 1) {
                double r_min = 1000.0; int k = 0;
                for (int j = 1; j <= nbonds[i]; ++j) {
                    int ii = ibonds[j][i];
                    double r = distance_(&ii, &i);
                    if (r < r_min && r > 0.95) { r_min = r; k = ibonds[j][i]; }
                }
                if (k > 0) {
                    for (int j = 1; j <= nbonds[i]; ++j) {
                        int l = ibonds[j][i]; int m = 1;
                        while (m <= nbonds[l] && ibonds[m][l] != i) m++;
                        nbonds[l]--;
                        for (int n = m; n <= nbonds[l]; ++n) ibonds[n][l] = ibonds[n+1][l];
                    }
                    nbonds[i] = 1; ibonds[1][i] = k;
                    nbonds[k]++; ibonds[nbonds[k]][k] = i;
                }
            }
        }
    }
}
