// ============================================================================
// 【中文导读】deri21.cpp —— 最小二乘正交基（对应 deri21.F90）
// ----------------------------------------------------------------------------
// 对 A'A 做打包矩阵本征分解，取累积方差达到阈值的前 ncut 个方向，构造
// B = A*V*diag(|lambda|^-1/2)（截断后正交化）；供 deri2 的 CI 梯度/Lagrange
// 步骤使用。布局：a/b 为 0-based 列主序（对应 F90 a(minear,nvar_nvo)），
// 与 C++ 内部 1-based 打包数组相互转换时要注意列主序。
// ============================================================================

// deri21.cpp — C++ translation of MOPAC 2016 "deri21.F90".
// Least-squares orthonormal basis: eigendecomposition of A'A (packed), then
// B = A * V * diag(|lambda|^-1/2) truncated at cumulative variance cutoff.
// Layout: a/b column-major 0-based (F90 a(minear,nvar_nvo), b(minear,ncut)).
#include "deri21.h"
#include <cmath>
#include <vector>
#include "mtxmc.h"
#include "rsp.h"
#include "mxm.h"

void deri21(const double* a, int nvar_nvo, int minear, int ifirst,
            double* vnert, double* pnert, double* b, int& ncut) {
    const double cutoff = 0.85;
    double sum2 = 0.0;
    // VNERT = A' * A, packed lower triangle (size nvar_nvo*(nvar_nvo+1)/2)
    std::vector<double> work((nvar_nvo * nvar_nvo + 1) * 2, 0.0);
    mtxmc(a, nvar_nvo, a, minear, work.data());
    int np = (nvar_nvo * (nvar_nvo + 1)) / 2;
    for (int k = 0; k < np; ++k) work[k] = -work[k];
    if (std::fabs(work[0]) < 1e-28 && nvar_nvo == 1) {
        pnert[0] = std::sqrt(-work[0]);
        work[0] = 1e15;
        vnert[0] = 1.0;
        ncut = 1;
    } else {
        rsp(work.data(), nvar_nvo, pnert, vnert);
        // eigenvalues ascending (LAPACK convention, matching MOPAC rsp)
        double sum = 0.0;
        for (int i = 0; i < nvar_nvo; ++i) sum -= pnert[i];
        int l = 0;
        int i = 0;
        for (; i < ifirst; ++i) {
            sum2 -= pnert[i] / sum;
            pnert[i] = std::sqrt(std::fabs(pnert[i]));
            for (int q = 0; q < nvar_nvo; ++q) work[l + q] = vnert[l + q] / pnert[i];
            l += nvar_nvo;
            if (sum2 < cutoff) continue;
            ncut = i + 1;
            break;
        }
        if (i >= ifirst) ncut = ifirst;
    }
    // B(minear,ncut) = A(minear,nvar_nvo) * work(nvar_nvo,ncut)
    mxm(a, minear, work.data(), nvar_nvo, b, ncut);
}
