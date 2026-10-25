// cnvg.cpp — C++ translation of MOPAC 2016 "cnvg.F90".
// ============================================================================
// 【中文导读】cnvg() —— SCF 密度收敛判据 + Camp-King 外推
// ----------------------------------------------------------------------------
// 输入：pnew 为本次迭代的新密度（打包），p 为上一次密度，p1 保存再上一次的
// 对角元；输出 pl = max|Δp|（SCF 收敛判定的核心指标）。
// 流程：
//   1) 遍历上三角：累计 pl，并按需计算二阶差分，得到外推系数 fac（隔轮
//      extrap）；把 p 更新为 pnew，p1 记录旧值；
//   2) 第二遍遍历：非对角元做 p += fac*(p-p1) 的外推，对角元限幅到
//      [0, rhfuhf]（UHF 时 1，否则 2），同时计算电子数归一化；
//   3) 若总电子数偏差 >1e-5，按比例缩放密度并重新限幅，循环至收敛。
// 约定：打包下三角 1-based；rhfuhf 由 numcal 首次调用时按 UHF 关键字设置。
// ============================================================================

#include "cnvg.h"

#include <algorithm>
#include <cmath>
#include <string>

#include "molkst_C.h"

using namespace molkst_C;

namespace {
double rhfuhf = 2.0;
int icalcn = -1;
}

void cnvg(std::vector<double>& pnew, std::vector<double>& p,
          std::vector<double>& p1, int niter, double& pl) {
    if (icalcn != numcal) {
        icalcn = numcal;
        rhfuhf = (keywrd.find("UHF") != std::string::npos) ? 1.0 : 2.0;
    }

    // pl = 最大密度变化量；damp 限制 >3 轮后的非对角变化，防止振荡。

    pl = 0.0;
    double faca = 0.0, damp = 1e10, facb = 0.0, fac = 0.0;
    if (niter > 3) damp = 0.05;
    bool extrap = (niter % 3) != 0;

    double sum1 = 0.0;
    int k = 0;
    for (int i = 1; i <= norbs; ++i) {
        k += i;
        double a = pnew[k];
        sum1 += a;
        double sa = std::abs(a - p[k]);
        if (sa > pl) pl = sa;
        if (!extrap) {
            faca += sa * sa;
            facb += (a - 2.0 * p[k] + p1[i]) * (a - 2.0 * p[k] + p1[i]);
        }
        p1[i] = p[k];
        p[k] = a;
    }

    // Camp-King 型外推系数：用相邻两轮的一阶/二阶差分比估计，满足条件才启用。

    if (facb > 1e-10 && faca < 100.0 * facb) fac = std::sqrt(faca / facb);

    int ie = 0;
    double sum2 = 0.0;
    for (int i = 1; i <= norbs; ++i) {
        int ii = i - 1;
        for (int j = 1; j <= ii; ++j) {
            ie++;
            double a = pnew[ie];
            p[ie] = a + fac * (a - p[ie]);
            pnew[ie] = p[ie];
        }
        ie++;
        if (std::abs(p[ie] - p1[i]) > damp) {
            p[ie] = p1[i] + std::copysign(damp, p[ie] - p1[i]);
        } else {
            p[ie] = p[ie] + fac * (p[ie] - p1[i]);
        }
        p[ie] = std::min(rhfuhf, std::max(p[ie], 0.0));
        sum2 += p[ie];
        pnew[ie] = p[ie];
    }

    // 电子数守恒循环：缩放密度后重新限幅，直到总电子数与 sum1 一致。

    double sum0 = sum1;
    while (true) {
        double sum = (sum2 > 1e-3) ? sum1 / sum2 : 0.0;
        sum1 = sum0;
        if (sum2 < 1e-3 || std::abs(sum - 1.0) < 1e-5) break;
        sum2 = 0.0;
        for (int i = 1; i <= norbs; ++i) {
            int j = i * (i + 1) / 2;
            p[j] = p[j] * sum + 1e-20;
            p[j] = std::max(p[j], 0.0);
            if (p[j] > rhfuhf) {
                p[j] = rhfuhf;
                sum1 -= rhfuhf;
            } else {
                sum2 += p[j];
            }
            pnew[j] = p[j];
        }
    }
}
