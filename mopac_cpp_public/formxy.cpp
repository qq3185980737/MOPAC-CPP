// ============================================================================
// 【中文导读】formxy.cpp —— CI xy 超向量组装（formxy.F90）
// ----------------------------------------------------------------------------
// 由 MO 系数 ca/cb 与权重 w 构造 CI 双电子超向量 xy（wca/wcb），并维护计数 kr；
// 1-based 打包下三角，供 dijkl/diagi 等 CI 能量与梯度使用。
// ============================================================================

// formxy.cpp — C++ translation of MOPAC 2016 "formxy.F90".

#include "formxy.h"

#include <vector>

void formxy(const std::vector<double>& w, int& kr,
            std::vector<double>& wca, std::vector<double>& wcb,
            const std::vector<double>& ca, const std::vector<double>& cb,
            int na, int nb) {
    std::vector<int> in(46, 0);
    in[1] = 1; in[10] = 4; in[45] = 9;
    int ij = 0, n1 = 0;
    int nna = in[na], nnb = in[nb];
    for (int i = 1; i <= nna; ++i) {
        double aa = 1.0;
        for (int j = 1; j <= i; ++j) {
            ++n1;
            if (i == j) aa = 0.5;
            ++ij;
            double sum = 0.0;
            int kl = 0, n2 = 0;
            for (int k = 1; k <= nnb; ++k) {
                double bb = 1.0;
                for (int l = 1; l <= k; ++l) {
                    ++n2;
                    if (k == l) bb = 0.5;
                    ++kl;
                    sum += cb[kl] * w[(n1 - 1) * nb + n2 - 1] * bb;
                }
            }
            wca[ij] += sum * aa;
        }
    }
    ij = 0; n1 = 0;
    for (int i = 1; i <= nnb; ++i) {
        double aa = 1.0;
        for (int j = 1; j <= i; ++j) {
            ++n1;
            if (i == j) aa = 0.5;
            ++ij;
            double sum = 0.0;
            int kl = 0, n2 = 0;
            for (int k = 1; k <= nna; ++k) {
                double bb = 1.0;
                for (int l = 1; l <= k; ++l) {
                    ++n2;
                    if (k == l) bb = 0.5;
                    ++kl;
                    sum += ca[kl] * w[(n2 - 1) * nb + n1 - 1] * bb;
                }
            }
            wcb[ij] += sum * aa;
        }
    }
    kr += na * nb;
}
