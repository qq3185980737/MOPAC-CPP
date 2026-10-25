// ============================================================================
// 【中文导读】dot.cpp —— 1-based 向量点积（dot.F90）
// ----------------------------------------------------------------------------
// dot(x,y,n) = sum_{i=1..n} x[i]*y[i]，x/y 为 1-based（下标 0 作填充）。
// 注意与 BLAS ddot 的区别：这里显式保留 1-based 约定，供 SCF/几何等复用。
// ============================================================================

// dot.cpp — C++ translation of MOPAC 2016 "dot.F90".
// Fortran x(n), y(n) are 1-based; convention keeps padding at [0].
#include "dot.h"

double dot(const std::vector<double>& x, const std::vector<double>& y, int n) {
    double s = 0.0;
    for (int i = 1; i <= n; ++i) s += x[i] * y[i];
    return s;
}
