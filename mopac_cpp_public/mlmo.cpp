// ============================================================================
// 【中文导读】mlmo.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mlmo.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mlmo.cpp — C++ translation of MOPAC 2016 "mlmo.F90".
// Bookkeeping routine for adding one LMO to the MOZYME occupied/virtual arrays.
// Fortran 1-based indexing preserved; pointer args are base + 0 (padding).

#include "mlmo.h"

#include <algorithm>

#include "molkst_C.h"
#include "MOZYME_C.h"

void mlmo(int& locc, int& lvir, int ii, int jj, int& nf_loc, int& ne,
          int& nocc, int& nvir, int* iz, int* ib, int* nce, int* ncf,
          int* ncocc, int* ncvir, const int* iorbs, int* icocc, int* icvir,
          double* cocc, double* cvir) {
    int nes = ne;
    int nfs = nf_loc;
    int iocc = locc;
    int ivir = lvir;

    if (ii != 0) {
        // OCCUPIED M.O.
        iz[ii] = iz[ii] - 1;
        if (jj == 0) {
            iz[ii] = iz[ii] - 1;
        }
        ib[ii] = ib[ii] - 1;
        nocc = nocc + 1;
        ncocc[nocc] = locc;
        locc = locc + iorbs[ii];
        nf_loc = nf_loc + 1;
        icocc[nf_loc] = ii;
        ncf[nocc] = 1;
    }
    if (jj != 0) {
        // VIRTUAL M.O.
        iz[jj] = iz[jj] - 1;
        if (ii == 0) {
            iz[jj] = iz[jj] + 1;
        }
        ib[jj] = ib[jj] - 1;
        nvir = nvir + 1;
        ncvir[nvir] = lvir;
        lvir = lvir + iorbs[jj];
        ne = ne + 1;
        nce[nvir] = 1;
        if (ii != 0) {
            icvir[ne] = ii;
            ncf[nocc] = 2;
            nce[nvir] = 2;
        } else {
            icvir[ne] = jj;
        }
    }
    if (ii != 0 && jj != 0) {
        // OCCUPIED AND VIRTUAL M.O.
        nf_loc = nf_loc + 1;
        icocc[nf_loc] = jj;
        ne = ne + 1;
        icvir[ne] = jj;
        locc = locc + iorbs[jj];
        lvir = lvir + iorbs[ii];
    }

    // Assign space for LMO expansion.
    int j = std::min(molkst_C::numat * 2, MOZYME_C::ipad2);
    int k = std::min(molkst_C::norbs * 2, MOZYME_C::ipad4);
    if (ii != 0) {
        nf_loc = nfs + j;
        for (int i = locc + 1; i <= iocc + k; ++i) {
            cocc[i] = 0.0;
        }
        locc = iocc + k;
    }
    if (jj != 0) {
        ne = nes + j;
        for (int i = lvir + 1; i <= ivir + k; ++i) {
            cvir[i] = 0.0;
        }
        lvir = ivir + k;
    }
}
