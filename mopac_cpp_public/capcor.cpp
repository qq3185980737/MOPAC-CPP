// capcor.cpp — C++ translation of MOPAC 2016 "capcor.F90".
// ============================================================================
// 【中文导读】capcor() —— 封端原子（CAP）能量修正
// ----------------------------------------------------------------------------
// 对应 capcor.F90。MOPAC 用原子序数 102 表示“封端伪原子”（CAP，用于截断
// 蛋白/片段时补足悬键）。本函数累加 CAP 原子相关的 p_ij*h_ij 项并乘 -2
// （闭壳层双占据约定），返回电子能量修正；iter() 在 SCF 收敛后若 capps
// 为真就把 capcor 加到电子能量上。
// 两个分支：CAP 原子自身的轨道块（ni==102），以及其它原子与 CAP 原子的
// 交叉块（nat[j]==102）。
// ============================================================================

#include "capcor.h"

#include "molkst_C.h"

// capcor(): 遍历原子累加 CAP(102) 相关的 p*h 项，返回 -2*sum。

double capcor(const std::vector<int>& nat, const std::vector<int>& nfirst,
              const std::vector<int>& nlast, const std::vector<double>& p,
              const std::vector<double>& h) {
    double sum = 0.0;
    for (int i = 1; i <= molkst_C::numat; ++i) {
        int ni = nat[i];
        int il = nfirst[i];
        int iu = nlast[i];
        if (ni == 102) {
            int j = (nlast[i] * (nlast[i] + 1)) / 2;
            int ii = iu - 1;
            for (int k = 1; k <= ii; ++k) {
                --j;
                sum += p[j] * h[j];
            }
        } else {
            for (int j = 1; j <= i; ++j) {
                int jl = nfirst[j];
                if (nat[j] != 102) continue;
                for (int k = il; k <= iu; ++k) {
                    int kk = (k * (k - 1)) / 2 + jl;
                    sum += p[kk] * h[kk];
                }
            }
        }
    }
    return -2.0 * sum;
}
