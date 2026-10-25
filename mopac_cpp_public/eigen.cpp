// ============================================================================
// 【中文导读】eigen.cpp —— 矩阵对角化/本征值处理
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 eigen.F90 的 C++ 移植（或同源子程序）。
// 本文件包含：
//   * eigen_limits() —— 解析 VECTORS 关键字；
//   * eigen()        —— 非 MOZYME 的 MO 向量输出（仍为占位，非 MOZYME 的 VEC
//                       输出由 writmo 现有代码处理）；
//   * FlmoFromFao() / lmo_to_eigenvectors() —— MOZYME 压缩 LMO 在 LMO 基下
//     构造 Fock 矩阵并对角化，得到"本征向量"（供 MULLIK/GRAPH 使用）。
// 约定：保持 Fortran 1-based 索引与列主序布局。
// ============================================================================

// eigen.cpp — C++ translation of MOPAC 2016 "eigen.F90".

#include "eigen.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "molkst_C.h"
#include "common_arrays_C.h"
#include "MOZYME_C.h"
#include "ijbo.h"
#include "rsp.h"

using namespace molkst_C;
using namespace MOZYME_C;
using common_arrays_C::f;
using common_arrays_C::nfirst;
using common_arrays_C::nlast;

void eigen_limits(int& print_nocc, int& print_nvir) {
    int noccupied = nelecs / 2;
    int nvirtual = norbs - noccupied;
    if (keywrd.find(" ALLVEC") != std::string::npos) {
        print_nocc = noccupied; print_nvir = nvirtual;
    } else if (keywrd.find(" VECTORS(") != std::string::npos ||
               keywrd.find(" VECTORS=(") != std::string::npos) {
        print_nocc = 8; print_nvir = 8;
    } else {
        print_nocc = 8; print_nvir = 8;
    }
    print_nocc = std::max(0, std::min(print_nocc, noccupied));
    print_nvir = std::max(0, std::min(print_nvir, nvirtual));
}

void eigen(bool, bool) {}

// ============================================================================
// FlmoFromFao(): F90 eigen.F90:192-379 —— 在 LMO 基下构造 Fock 矩阵。
// ============================================================================
void FlmoFromFao(int nlmo, const std::vector<int>& ncx, const std::vector<int>& nncx,
                 const std::vector<int>& ncxxx, int /*n_dim*/,
                 std::vector<char>& latoms, std::vector<double>& ws,
                 const std::vector<int>& icxxx, int /*icxxx_dim*/,
                 const std::vector<double>& cxxx, int /*cxxx_dim*/,
                 std::vector<double>& ac, std::vector<double>& fmo_expanded) {
    const double flim = 1.0e-8;
    const double cutoff = 1.0e-8;
    std::fill(fmo_expanded.begin(), fmo_expanded.end(), 0.0);

    int k_fmo = 1;
    for (int i = 1; i <= nlmo; ++i) {
        const int loopi = ncxxx[i];

        // Weight on each atom from this LMO.
        for (int j = 1; j <= numat; ++j) latoms[j] = 0;
        int l = 0;
        for (int j = nncx[i] + 1; j <= nncx[i] + ncx[i]; ++j) {
            const int j1 = icxxx[j];
            double sum = 0.0;
            for (int k = l + 1; k <= l + iorbs[j1]; ++k)
                sum += cxxx[k + loopi] * cxxx[k + loopi];
            ac[j1] = sum;
            latoms[j1] = 1;
            l += iorbs[j1];
        }

        // Evaluate F * C.
        for (int jj = nncx[i] + 1; jj <= nncx[i] + ncx[i]; ++jj) {
            const int j1 = icxxx[jj];
            for (int jx = 1; jx <= iorbs[j1]; ++jx)
                ws[nfirst[j1] + jx - 1] = 0.0;
            int kl = loopi;
            for (int kk = nncx[i] + 1; kk <= nncx[i] + ncx[i]; ++kk) {
                const int k1 = icxxx[kk];
                const int kj = lijbo ? nijbo[k1][j1] : ijbo(k1, j1);
                if (kj >= 0 && ac[k1] > cutoff) {
                    if (iorbs[k1] == 1 && iorbs[j1] == 1) {
                        ws[nfirst[j1]] += f[kj + 1] * cxxx[kl + 1];
                    } else if (k1 > j1) {
                        int ii = kj;
                        for (int i4 = 1; i4 <= iorbs[k1]; ++i4)
                            for (int jx = 1; jx <= iorbs[j1]; ++jx) {
                                ++ii;
                                ws[nfirst[j1] + jx - 1] += f[ii] * cxxx[kl + i4];
                            }
                    } else if (k1 < j1) {
                        int ii = kj;
                        for (int jx = 1; jx <= iorbs[j1]; ++jx)
                            for (int i4 = 1; i4 <= iorbs[k1]; ++i4) {
                                ++ii;
                                ws[nfirst[j1] + jx - 1] += f[ii] * cxxx[kl + i4];
                            }
                    } else {
                        for (int jx = 1; jx <= iorbs[j1]; ++jx)
                            for (int i4 = 1; i4 <= iorbs[j1]; ++i4) {
                                int ii = kj + (jx * (jx - 1)) / 2 + i4;
                                if (i4 > jx) ii = kj + (i4 * (i4 - 1)) / 2 + jx;
                                ws[nfirst[j1] + jx - 1] += f[ii] * cxxx[kl + i4];
                            }
                    }
                }
                kl += iorbs[k1];
            }
        }

        // Evaluate C * (F * C).
        const int i1 = icxxx[nncx[i] + 1];
        const int i2 = (ncx[i] > 1) ? icxxx[nncx[i] + 2] : i1;
        for (int j = 1; j <= i; ++j) {
            const int j1 = icxxx[nncx[j] + 1];
            const int j2 = (ncx[j] > 1) ? icxxx[nncx[j] + 2] : j1;
            int i1j1, i1j2, i2j1, i2j2;
            if (lijbo) {
                i1j1 = nijbo[i1][j1]; i1j2 = nijbo[i1][j2];
                i2j1 = nijbo[i2][j1]; i2j2 = nijbo[i2][j2];
            } else {
                i1j1 = ijbo(i1, j1); i1j2 = ijbo(i1, j2);
                i2j1 = ijbo(i2, j1); i2j2 = ijbo(i2, j2);
            }
            if (i1j1 >= 0 || i1j2 >= 0 || i2j1 >= 0 || i2j2 >= 0) {
                double sum = 0.0;
                if (i1j1 >= 0) sum = std::fabs(f[i1j1 + 1]);
                if (i1j2 >= 0) sum += std::fabs(f[i1j2 + 1]);
                if (i2j1 >= 0) sum += std::fabs(f[i2j1 + 1]);
                if (i2j2 >= 0) sum += std::fabs(f[i2j2 + 1]);
                if (sum >= flim) {
                    const int loopj = ncxxx[j];
                    double sum2 = 0.0;
                    int kl = 0;
                    for (int kk = nncx[j] + 1; kk <= nncx[j] + ncx[j]; ++kk) {
                        const int k1 = icxxx[kk];
                        if (!latoms[k1]) {
                            kl += iorbs[k1];
                        } else {
                            for (int k = nfirst[k1]; k <= nlast[k1]; ++k) {
                                ++kl;
                                sum2 += ws[k] * cxxx[kl + loopj];
                            }
                        }
                    }
                    fmo_expanded[k_fmo] = sum2;
                }
            }
            ++k_fmo;
        }
    }
}

// ============================================================================
// lmo_to_eigenvectors(): F90 eigen.F90:438-499 —— 把 LMO 基 Fock 对角化，
// 再展开回 AO 基，得到常规本征向量（供 Mulliken 分析用）。
// eig_offset / col_offset 支持虚拟集从 noccupied+1 开始。
// ============================================================================
void lmo_to_eigenvectors(int nlmo, const std::vector<int>& ncx, const std::vector<int>& nncx,
                         const std::vector<int>& ncxxx, int n_dim,
                         const std::vector<int>& icxxx, int icxxx_dim,
                         const std::vector<double>& cxxx, int cxxx_dim,
                         std::vector<double>& eigs_out, int eig_offset,
                         std::vector<std::vector<double>>& c, int col_offset) {
    std::vector<char> latoms(numat + 1, 0);
    std::vector<double> ws(norbs + 1, 0.0);
    std::vector<double> ac(numat + 1, 0.0);
    std::vector<double> fmo_expanded((size_t)nlmo * nlmo + 1, 0.0);

    FlmoFromFao(nlmo, ncx, nncx, ncxxx, n_dim, latoms, ws, icxxx, icxxx_dim,
                cxxx, cxxx_dim, ac, fmo_expanded);

    // Diagonalise Fock matrix in the LMO basis.
    std::vector<double> vec_lmo((size_t)nlmo * nlmo + 1, 0.0);
    std::vector<double> eig_tmp(nlmo + 1, 0.0);
    rsp(&fmo_expanded[1], nlmo, &eig_tmp[1], &vec_lmo[1]);
    for (int i = 1; i <= nlmo; ++i) eigs_out[eig_offset + i] = eig_tmp[i];

    // Reconstruct eigenvectors in AO basis.
    for (int k = 1; k <= norbs; ++k)
        for (int j = 1; j <= nlmo; ++j) c[k][col_offset + j] = 0.0;

    std::vector<double> c_lmo(norbs + 1, 0.0);
    for (int i = 1; i <= nlmo; ++i) {
        std::fill(c_lmo.begin(), c_lmo.end(), 0.0);
        int ka = ncxxx[i];
        for (int jj = nncx[i] + 1; jj <= nncx[i] + ncx[i]; ++jj) {
            const int j = icxxx[jj];
            for (int k = nfirst[j]; k <= nlast[j]; ++k) {
                ++ka;
                c_lmo[k] = cxxx[ka];
            }
        }
        for (int k = 1; k <= norbs; ++k) {
            if (std::fabs(c_lmo[k]) <= 1.0e-10) continue;
            int l = i;
            for (int j = 1; j <= nlmo; ++j) {
                c[k][col_offset + j] += vec_lmo[l] * c_lmo[k];
                l += nlmo;
            }
        }
    }
}
