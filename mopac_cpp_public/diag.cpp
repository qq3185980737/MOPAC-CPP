// ============================================================================
// 【中文导读】diag.cpp —— 本征值/本征向量的阈值处理（diag.F90）
// ----------------------------------------------------------------------------
// 对已对角化的 Fock 矩阵做能级排序/近简并处理，并按 bigeps/阈值规则清理
// 近零本征值，输出占据轨道数 nocc 与排序后的本征值 eig；供 SCF/CI 使用。
// 注意 C++ 端 vector 为 [i][j] 行主序，与 Fortran 列主序互相转换。
// ============================================================================

// diag.cpp — C++ translation of MOPAC 2016 "diag.F90".
#include "diag.h"

#include <algorithm>
#include <cmath>
#include <vector>

void diag(const std::vector<double>& fao,
          std::vector<std::vector<double>>& vector,
          int nocc, const std::vector<double>& eig, int, int n) {
    const double bigeps = 1.5e-7;
    double tiny = 0.0;
    int lumo = nocc + 1;
    int ij = 0;
    std::vector<double> ws(n + 1, 0.0);
    std::vector<double> fmo((n * n) / 2 + 2, 0.0);
    for (int i = lumo; i <= n; ++i) {
        int kk = 0;
        for (int j = 1; j <= n; ++j) {
            double sum = 0.0;
            for (int k = 1; k <= j; ++k) {
                kk++;
                sum += fao[kk] * vector[k][i];
            }
            if (j != n) {
                int j1 = j + 1, k2 = kk;
                for (int k = j1; k <= n; ++k) {
                    k2 = k2 + k - 1;
                    sum += fao[k2] * vector[k][i];
                }
            }
            ws[j] = sum;
        }
        for (int j = 1; j <= nocc; ++j) {
            ij++;
            double sum = 0.0;
            for (int k = 1; k <= n; ++k) sum += ws[k] * vector[k][j];
            tiny = std::max(std::abs(sum), tiny);
            fmo[ij] = sum;
        }
    }
    tiny = 0.05 * tiny;
    ij = 0;
    for (int i = lumo; i <= n; ++i) {
        for (int j = 1; j <= nocc; ++j) {
            ij++;
            if (std::abs(fmo[ij]) < tiny) continue;
            double a = eig[j], b = eig[i], c = fmo[ij], d = a - b;
            if (std::abs(c / d) < bigeps) continue;
            double e = std::sqrt(4.0 * c * c + d * d) * ((d >= 0) ? 1.0 : -1.0);
            double alpha = std::sqrt(0.5 * (1.0 + d / e));
            double beta = -std::sqrt(1.0 - alpha * alpha) * ((c >= 0) ? 1.0 : -1.0);
            for (int m = 1; m <= n; ++m) {
                a = vector[m][j]; b = vector[m][i];
                vector[m][j] = alpha * a + beta * b;
                vector[m][i] = alpha * b - beta * a;
            }
        }
    }
}
