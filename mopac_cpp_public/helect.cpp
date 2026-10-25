#include <vector>
// helect.cpp
#include "helect.h"

// ============================================================================
// 【中文导读】helect() —— 计算电子能量（E = Σ_ij P_ij (H_ij + F_ij)）
// ----------------------------------------------------------------------------
// 对应 helect.F90。输入 P/H/F 均为 1-based 打包下三角（索引 j*(j-1)/2+i）。
// 实现用双重循环只遍历上三角：先累加对角元 ed，再累加严格上三角 ee，
// 最后返回 ee + ed/2（对角元只计一次，非对角元被对称地计两次再减半）。
// 注意：k 从 0 开始但第一次 ++k 后即为 Fortran 的 1，所以所有数组访问都按
// Fortran 下标；helect 返回 eV（MOPAC 内部能量单位）。
// ============================================================================

// Fortran packed lower-triangle arrays are 1-based; p/h/f[0] is padding.
double helect(int n, const std::vector<double>& p, const std::vector<double>& h, const std::vector<double>& f) { return helect(n, p.data(), h.data(), f.data()); }

// 核心实现（指针版本；vector 版本在上一行做 .data() 转发）。

double helect(int n, const double* p, const double* h, const double* f) {
    double ed = 0.0, ee = 0.0;
    int k = 0;  // Fortran k starts at 1 after first increment
    int nn = n + 1;
    for (int i = 2; i <= nn; ++i) {
        ++k;  // now Fortran k = 1, 2, ...
        int jj = i - 1;
        ed += p[k] * (h[k] + f[k]);
        if (i == nn) continue;
        if (jj > 0) {
            double s = 0.0;
            for (int t = 1; t <= jj; ++t) s += p[k+t] * (h[k+t] + f[k+t]);
            ee += s;
            k = jj + k;
        }
    }
    ee += 0.5 * ed;
    return ee;
}
