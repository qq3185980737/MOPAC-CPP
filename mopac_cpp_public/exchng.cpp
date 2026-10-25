// ============================================================================
// 【中文导读】exchng.cpp —— 标量/向量片段交换（exchng.F90）
// ----------------------------------------------------------------------------
// 交换两组标量 (a,c) <-> (b,d) 并按参数 t 更新 q；同时交换 x/y 向量片段
// （1-based，长度 n）。用于 MO 排序/对称化中的行列交换。
// ============================================================================

// exchng.cpp — C++ translation of MOPAC 2016 "exchng.F90".

#include "exchng.h"

void exchng(double a, double& b, double c, double& d,
            double t, double& q, const std::vector<double>& x,
            std::vector<double>& y, int n) {
    b = a;
    d = c;
    q = t;
    for (int i = 1; i <= n; ++i) y[i] = x[i];
}
