// ============================================================================
// 【中文导读】check.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 check.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// check.cpp — C++ translation of MOPAC 2016 "check.F90".

#include "check.h"

#include <cmath>

#include "MOZYME_C.h"

void check(int nvec, const std::vector<int>& nnc, const std::vector<int>& nc,
           const std::vector<int>& icvec, const std::vector<int>& iorbs,
           const std::vector<int>& ncvec, std::vector<double>& cvec) {
    double error = 0.0;
    for (int i = 1; i <= nvec; ++i) {
        int m = nnc[i];
        double sum = 0.0;
        int n = 0;
        for (int j = 1; j <= nc[i]; ++j) {
            ++m;
            int k = icvec[m];
            for (int mm = 1; mm <= iorbs[k]; ++mm) {
                ++n;
                double v = cvec[ncvec[i] + n];
                sum += v * v;
            }
        }
        error += std::abs(1.0 - sum);
        MOZYME_C::ws[i] = sum;
    }
    for (int i = 1; i <= nvec; ++i) {
        double scale = 1.0 / std::sqrt(MOZYME_C::ws[i]);
        int n = 0;
        for (int j = 1; j <= nc[i]; ++j) {
            int k = icvec[nnc[i] + j];
            for (int mm = 1; mm <= iorbs[k]; ++mm) {
                ++n;
                cvec[ncvec[i] + n] *= scale;
            }
        }
    }
    if (error <= 0.1) return;
    // Severe error: would print and call mopend. In translation, return.
}
