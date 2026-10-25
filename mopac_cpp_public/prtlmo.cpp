// ============================================================================
// 【中文导读】prtlmo.cpp —— MOZYME 局域分子轨道（LMO）打印
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 prtlmo.F90（prtlmo + prtlmn）。按 LMO 能量排序后，把每个
// LMO 在原子上的 |psi|^2（乘 10000 取整）按贡献从大到小列出，最多 20 个原子。
// 数据来自 MOZYME_C 的 cocc/cvir/icocc/icvir/ncocc/ncvir/nncf/nnce 及共同
// 数组 eigs/nat；保持 Fortran 1-based 索引。
// ============================================================================

// prtlmo.cpp
#include "prtlmo.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "elemts_C.h"
#include "eigen.h"
using namespace MOZYME_C;
using common_arrays_C::nat;
using common_arrays_C::eigs;
using molkst_C::norbs;
using MOZYME_C::iorbs;
using MOZYME_C::noccupied;
using MOZYME_C::nvirtual;
using elemts_C::elemnt;

// prtlmn: Fortran prtlmo.F90:40-160.
// isort/eigs 用偏移量传递：占据集从 1 开始，虚拟集从 noccupied+1 开始。
static void prtlmn(const std::vector<int>& nncx, const std::vector<int>& icxxx,
                   const std::vector<int>& ncxxx, const std::vector<double>& cxxx,
                   const std::vector<int>& ncx, std::vector<int>& isort, int isort_off,
                   const std::vector<double>& eigv, int eig_off,
                   int mmos, int i_offset) {
    std::vector<double> eigs_temp(mmos + 1, 0.0);
    for (int i = 1; i <= mmos; ++i) eigs_temp[i] = eigv[eig_off + i];

    // Construct map of LMO energy levels.
    for (int i = 1; i <= mmos; ++i) {
        double eig_min = 1.0e7;
        int ii = 1;
        for (int j = 1; j <= mmos; ++j)
            if (eigs_temp[j] < eig_min) { eig_min = eigs_temp[j]; ii = j; }
        isort[isort_off + i] = ii;
        eigs_temp[ii] = 1.0e8;
    }

    // Put LMO contributions per atom into an atom-based array; JAT holds the
    // number of important atoms in each LMO.
    const int cap = std::max(1, noccupied * 20);
    std::vector<double> w(cap + 1, 0.0);
    std::vector<int> iscrch(cap + 1, 0);
    std::vector<int> jat(norbs + 1, 0);
    const double cnst = 2.0e-4;
    int j = 0;
    for (int iunsrt = 1; iunsrt <= mmos; ++iunsrt) {
        int i = isort[isort_off + iunsrt];
        int n = ncxxx[i];
        int k = 0;
        for (int ii = nncx[i] + 1; ii <= nncx[i] + ncx[i]; ++ii) {
            int l = icxxx[ii];
            int nj = iorbs[l];
            double sum = 0.0;
            for (int m = 1; m <= nj; ++m) {
                const double v = cxxx[m + n];
                sum += v * v;
            }
            if (sum > cnst) {
                ++j;
                ++k;
                iscrch[j] = l;
                w[j] = sum * 2.0;
            }
            n += nj;
        }
        jat[i] = k;
    }

    int ju = 0;
    std::vector<double> xbig(101, 0.0);
    std::vector<int> ibig(101, 0), jbig(101, 0);
    for (int iunsrt = 1; iunsrt <= mmos; ++iunsrt) {
        int i = isort[isort_off + iunsrt];
        int jl = ju + 1;
        ju = jl + jat[i] - 1;
        int k = 0;
        for (int jj = jl; jj <= ju; ++jj) {
            ++k;
            xbig[k] = w[jj];
            ibig[k] = iscrch[jj];
        }
        int kk = std::min(20, k);
        int loop = 1;
        for (loop = 1; loop <= kk; ++loop) {
            double sum = -1.0;
            int l = 0, lj = 0;
            for (int jj = loop; jj <= k; ++jj)
                if (xbig[jj] > sum) { sum = xbig[jj]; l = ibig[jj]; lj = jj; }
            if (sum < cnst) break;
            xbig[lj] = xbig[loop];
            ibig[lj] = ibig[loop];
            jbig[loop] = (int)(sum / cnst);
            ibig[loop] = l;
        }
        std::printf("%5d%9.4f", iunsrt + i_offset, eigv[eig_off + i]);
        for (int jj = 1; jj < loop; ++jj) {
            std::string sym = elemnt[nat[ibig[jj]]];
            while (!sym.empty() && sym.back() == ' ') sym.pop_back();
            std::printf("%5d %2s%5d", ibig[jj], sym.c_str(), jbig[jj]);
        }
        std::printf("\n");
    }
    (void)j;
}

// ----------------------------------------------------------------------------
// LMO 能量：官方 prtlmo 打印的 `eigs` 是 MOZYME SCF 中 DIAGG1 写入的
// LMO 期望值 <i|F|i>（不是 canonical 本征值）。这里用 FlmoFromFao 在 LMO 基
// 下构造 Fock 并取对角元，避免把整条 diagg1 路径再搬一遍。
// ----------------------------------------------------------------------------
static void compute_lmo_energies(int nlmo, const std::vector<int>& ncx,
                                 const std::vector<int>& nncx,
                                 const std::vector<int>& ncxxx,
                                 const std::vector<int>& icxxx,
                                 const std::vector<double>& cxxx,
                                 std::vector<double>& out) {
    out.assign((size_t)nlmo + 1, 0.0);
    if (nlmo <= 0) return;
    std::vector<char> latoms((size_t)molkst_C::numat + 1, 0);
    std::vector<double> ws((size_t)molkst_C::norbs + 1, 0.0);
    std::vector<double> ac((size_t)molkst_C::numat + 1, 0.0);
    std::vector<double> fmo((size_t)nlmo * nlmo + 1, 0.0);
    FlmoFromFao(nlmo, ncx, nncx, ncxxx, nlmo, latoms, ws, icxxx,
                (int)icxxx.size() - 1, cxxx, (int)cxxx.size() - 1, ac, fmo);
    for (int i = 1; i <= nlmo; ++i) out[i] = fmo[(size_t)i * (i + 1) / 2];
}

void prtlmo() {
    if (isort.size() < (size_t)norbs + 1) isort.assign((size_t)norbs + 1, 0);
    std::printf("  LOCALIZED MOLECULAR ORBITALS\n");
    std::printf("\n (PSI**2 IS 10000 TIMES THE SQUARE OF THE WAVE-FUNCTION ON AN ATOM)\n\n");

    // Occupied/virtual LMO energies: official prtlmo prints eigs(1:noccupied)
    // and eigs(noccupied+1:) as left by the SCF diagonalisation (the values
    // MOZYME keeps in eigs), not a fresh Fock expectation value.
    std::vector<double> eocc((size_t)noccupied + 1, 0.0), evir((size_t)nvirtual + 1, 0.0);
    for (int i = 1; i <= noccupied; ++i) eocc[i] = eigs[i];
    for (int i = 1; i <= nvirtual; ++i) evir[i] = eigs[noccupied + i];

    // Print the occupied set of LMOs.
    std::printf("\n%10s Occupied Set\n\n", "");
    std::printf("  LMO   ENERGY    ATOM PSI**2  ATOM PSI**2  ATOM PSI**2  ATOM PSI**2  ATOM PSI**2\n");
    prtlmn(nncf, icocc, ncocc, cocc, ncf, isort, 0, eocc, 0, noccupied, 0);

    // Print the virtual set of LMOs.
    std::printf("\n%10s Virtual Set\n\n", "");
    std::printf("  LMO   ENERGY    ATOM PSI**2  ATOM PSI**2  ATOM PSI**2  ATOM PSI**2  ATOM PSI**2\n");
    prtlmn(nnce, icvir, ncvir, cvir, nce, isort, noccupied, evir, 0,
           nvirtual, noccupied);
}
