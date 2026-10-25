// ============================================================================
// 【中文导读】mamult.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mamult.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mamult.cpp
#include "mamult.h"
#ifdef _OPENMP
#include <omp.h>
#endif
void mamult(const double* a, const double* b, double* c, int n, double one) {
    // Each packed-triangle output element c[l] is written by exactly one
    // (i,j) pair and its inner k-summation order is unchanged, so the outer
    // (i,j) loop is embarrassingly parallel and bit-level identical to the
    // serial version.  Parallel only when OMP threads > 1 (bit-level default).
#ifdef _OPENMP
    if (omp_get_max_threads() > 1 && n > 200) {
#pragma omp parallel for schedule(static)
        for (int i = 1; i <= n; ++i) {
            int ii = ((i - 1) * i) / 2;
            for (int j = 1; j <= i; ++j) {
                int jj = ((j - 1) * j) / 2;
                int lidx = ii + j - 1;  // 0-based packed index (serial l-1)
                double sum = 0.0;
                for (int k = 1; k <= j; ++k) sum += a[ii + k] * b[jj + k];
                for (int k = j + 1; k <= i; ++k) sum += a[ii + k] * b[((k - 1) * k) / 2 + j];
                for (int k = i + 1; k <= n; ++k) {
                    int kk = (k * (k - 1)) / 2;
                    sum += a[kk + i] * b[kk + j];
                }
                c[lidx] = sum + one * c[lidx];
            }
        }
        return;
    }
#endif
    int l = 0;
    for (int i = 1; i <= n; ++i) {
        int ii = ((i - 1) * i) / 2;
        for (int j = 1; j <= i; ++j) {
            int jj = ((j - 1) * j) / 2;
            ++l;
            double sum = 0.0;
            for (int k = 1; k <= j; ++k) sum += a[ii + k] * b[jj + k];
            for (int k = j + 1; k <= i; ++k) sum += a[ii + k] * b[((k - 1) * k) / 2 + j];
            for (int k = i + 1; k <= n; ++k) {
                int kk = (k * (k - 1)) / 2;
                sum += a[kk + i] * b[kk + j];
            }
            c[l] = sum + one * c[l];
        }
    }
}
