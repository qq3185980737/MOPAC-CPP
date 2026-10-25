// ============================================================================
// 【中文导读】dofs.cpp —— 自由能/态密度直方图（dofs.F90）
// ----------------------------------------------------------------------------
// 把能量参考点 eref 累积成直方图 dd（按 bottom..top 区间、m 个 bin），
// 用于热力学配分函数/自由能计算。1-based 循环，bin 宽度 range 由边界决定。
// ============================================================================

// dofs.cpp — C++ translation of MOPAC 2016 "dofs.F90".

#include "dofs.h"

#include <cstdio>

void dofs(std::vector<std::vector<double>>& eref, int mono3, int n,
          std::vector<double>& dd, int m, double bottom, double top) {
    for (int k = 0; k < (int)dd.size(); ++k) dd[k] = 0.0;
    double range = m / (top - bottom);
    for (int j = 1; j <= mono3; ++j)
        for (int i = 1; i <= n; ++i) {
            double x = eref[j][i];
            if (x < bottom || x > top) x = -1e7;
            eref[j][i] = (x - bottom) * range;
        }
    for (int ii = 1; ii <= mono3; ++ii)
        for (int i = 2; i <= n; ++i) {
            double b = eref[ii][i-1];
            if (b < 1) continue;
            double a = eref[ii][i];
            if (a < 1) continue;
            if (b > a) { double t=b; b=a; a=t; }
            int j = (int)b, k = (int)a;
            if (j == k) { dd[k] += 1.0; }
            else {
                double spread = 1.0 / (a - b + 1e-12);
                dd[j] += (j + 1 - b) * spread;
                dd[k] += (a - k) * spread;
                if (k != j + 1) {
                    for (int q = j+1; q <= k-1; ++q) dd[q] += spread;
                }
            }
        }
    double x = m / ((n - 1) * (top - bottom));
    for (int k = 0; k < (int)dd.size(); ++k) dd[k] *= x;
    std::printf(" NORMALIZED DENSITY OF STATES\n");
    range = m / (top - bottom);
    for (int i = 1; i <= m; ++i)
        std::printf("%9.2f%12.6f\n", bottom + (i - 0.5) / range, dd[i]);
}
