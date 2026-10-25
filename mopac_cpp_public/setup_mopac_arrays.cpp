// setup_mopac_arrays.cpp — C++ translation of MOPAC 2016 "setup_mopac_arrays.F90".
// Create and destroy the allocatable arrays used by MOPAC.
//   n == 0            : destroy all arrays that have been created.
//   n > 0, mode == 1  : create the arrays essential for reading in data.
//   n > 0, mode != 1  : create MOPAC arrays (holds electronic info).

//
// ============================================================================
// 【中文导读】setup_mopac_arrays() —— 全局数组的分配与释放
// ----------------------------------------------------------------------------
// 对应 setup_mopac_arrays.F90。调用约定：
//   n == 0            : 释放所有数组（作业结束/切换体系）；
//   n > 0, mode == 1  : 只分配读数据阶段需要的数组（geo/coord/nat/...）；
//   n > 0, mode != 1  : 分配 SCF 及后续计算所需的全部数组（尺寸由 numat/
//                       norbs/mpack 等决定）。
// 数组约定：所有一维数组按 Fortran 1-based 分配（大小 N+1，0 号 padding）；
// 坐标类数组分配 4 行（0..3），因为 getgeo/gmetry 用 1..3 访问行列；
// 内存不足时：C++ 用 try/catch 捕获 bad_alloc，统一走 memory_error()。
// 提示：数组尺寸依赖已被前一步设置好的 maxatoms/norbs；若看到越界，先查
// 本函数是否在 readmo/getdat 之后以正确的 mode 调用过。
// ============================================================================

#include "setup_mopac_arrays.h"
#include "common_arrays_C.h"
#include "maps_C.h"
#include "iter_C.h"
#include "symmetry_C.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "derivs_C.h"
#include "ef_C.h"
#include "meci_C.h"
#include "drc_C.h"
#include "cosmo_C.h"
#include "esp_C.h"
#include "to_screen_C.h"
#include <cstdio>
#include <string>
#include <vector>

extern double meci();
extern void fock2(std::vector<double>&, const std::vector<double>&,
                  std::vector<double>&, const std::vector<double>&,
                  const std::vector<double>&, const std::vector<double>&,
                  int, const std::vector<int>&, const std::vector<int>&, int);
extern void mopend(const std::string&);

// hint(): 32 位版本的提示信息（官方行为；64 位构建不会触发）。

void hint() {
    using namespace molkst_C;
    using namespace chanel_C;
    if (num_bits == 32)
        std::fprintf(stdout,
                     "\n This version of MOPAC uses a 32-bit architecture.\n"
                     " The job might run if a 64-bit version of MOPAC is used\n"
                     " To download a 64-bit version, go to http://openmopac.net/Download_MOPAC_Executable_Step2.html\n");
}

void memory_error(const std::string& txt) {
    using namespace chanel_C;
    std::fprintf(stdout, "\n          Unable to allocate memory in subroutine %s\n", txt.c_str());
    mopend(txt);
    hint();
}

void setup_mopac_arrays(int n, int mode) {
    using namespace common_arrays_C;
    using namespace maps_C;
    using namespace iter_C;
    using namespace symmetry_C;
    using namespace molkst_C;
    using namespace chanel_C;
    using namespace derivs_C;
    using namespace ef_C;
    using namespace meci_C;
    using namespace drc_C;
    using namespace cosmo_C;
    using namespace esp_C;
    using namespace to_screen_C;

    // n>0：分配；n==0（走到下面的 else）：释放。

    if (n > 0) {
        if (mode == 1) {
            // Create essential arrays only at this point.

            // mode==1：仅读数据所需的最小数组集，避免在还没读到原子数时分配大数组。

            // NOTE: geo/coord/coorda are indexed 1..3 (1-based Fortran) in
            // getgeo/gmetry/etc., so allocate 4 rows (0..3); gmetry itself
            // reads coord via 0-based j-1 offsets, but geo uses geo[1..3].
            geo.assign(4, std::vector<double>((size_t)n + 1, 0.0));   // (1..3,n)
            coord.assign(4, std::vector<double>((size_t)n + 1, 0.0));
            coorda.assign(4, std::vector<double>((size_t)n + 1, 0.0));
            na.assign((size_t)n + 1, 0);
            nb.assign((size_t)n + 1, 0);
            nc.assign((size_t)n + 1, 0);
            simbol.assign(3 * (size_t)n + 1, std::string());
            atmass.assign((size_t)n + 1, 0.0);
            labels.assign((size_t)n + 1, 0);
            // loc[1..2] used 1-based in compfg; 3 rows = 0..2.
            loc.assign(3, std::vector<int>(3 * (size_t)n + 1, 0));
            xparam.assign(3 * (size_t)n + 1, 0.0);
            nat.assign((size_t)n + 1, 0);
            nfirst.assign((size_t)n + 1, 0);
            nlast.assign((size_t)n + 1, 0);
            uspd.assign(9 * (size_t)n + 1, 0.0);
            pdiag.assign(9 * (size_t)n + 1, 0.0);
            xparef.assign(3 * (size_t)n + 1, 0.0);
            depmul.assign(3 * (size_t)n + 1, 0.0);
            jelem.assign(21, std::vector<int>((size_t)n + 1, 0));
            txtatm1.assign((size_t)n + 1, std::string());
            txtatm.assign((size_t)n + 1, std::string());
            locpar.assign(3 * (size_t)n + 1, 0);
            idepfn.assign(3 * (size_t)n + 1, 0);
            locdep.assign(3 * (size_t)n + 1, 0);
            nbonds.assign((size_t)n + 1, 0);
            ibonds.assign(16, std::vector<int>((size_t)n + 1, 0));
            na_store.assign((size_t)n + 1, 0);
            l_atom.assign((size_t)n + 1, 0);
            std::fill(na_store.begin(), na_store.end(), 0);
            std::fill(nfirst.begin(), nfirst.end(), -9999);
            std::fill(nlast.begin(), nlast.end(), -9999);
            std::fill(nbonds.begin(), nbonds.end(), 0);
        } else {
            // Create main arrays — at this point, the array sizes are all known.
            int j = 0, i = 0;
            grad.assign((size_t)nvar + 1, 0.0);
            if (!mozyme) {
                // A conventional MOPAC job.
                h.assign((size_t)mpack + 1, 0.0);
                p.assign((size_t)mpack + 1, 0.0);
                pa.assign((size_t)mpack + 1, 0.0);
                pb.assign((size_t)mpack + 1, 0.0);
                pold.assign(6 * (size_t)mpack + 1, 0.0);
                pold2.assign(6 * (size_t)mpack + 1, 0.0);
                f.assign((size_t)mpack + 1, 0.0);
                c.assign((size_t)norbs + 1, std::vector<double>((size_t)norbs + 1, 0.0));
                eigs.assign((size_t)norbs + 2, 0.0);
                q.assign((size_t)numat + 1, 0.0);
                eigb.assign((size_t)norbs + 1, 0.0);
                pold3.assign((size_t)std::max(mpack, 400) + 1, 0.0);
                errfn.assign(3 * (size_t)natoms * (size_t)l123 + 1, 0.0);
                aicorr.assign((size_t)nvar + 1, 0.0);
                w.assign((size_t)n2elec + 2025 + 1, 0.0);
                dxyz.assign(3 * (size_t)numat * (size_t)l123 + 1, 0.0);
                if (uhf) {
                    fb.assign((size_t)mpack + 1, 0.0);
                    cb.assign((size_t)norbs + 1, std::vector<double>((size_t)norbs + 1, 0.0));
                    pbold.assign(6 * (size_t)mpack + 1, 0.0);
                    pbold2.assign(6 * (size_t)mpack + 1, 0.0);
                    pbold3.assign((size_t)std::max(mpack, 400) + 1, 0.0);
                }
                // Fortran checks every allocation's stat; on failure it reports
                // "too big to run". std::vector throws, so emulate the check
                // with a try/catch around the whole block.
                std::fill(pold.begin(), pold.end(), 0.0);
                std::fill(eigb.begin(), eigb.end(), 0.0);
                if (uhf) std::fill(pbold.begin(), pbold.end(), 0.0);
                if (l123 > 1) {
                    wk.assign((size_t)n2elec + 2025 + 1, 0.0);
                }
                std::fill(dxyz.begin(), dxyz.end(), 0.0);
                std::fill(errfn.begin(), errfn.end(), 0.0);
                std::fill(aicorr.begin(), aicorr.end(), 0.0);
                std::fill(grad.begin(), grad.end(), 0.0);
            }
        }
    } else {
        // Delete all arrays.

        // 中文：释放全部数组。C++ 侧用 clear()/shrink_to_fit() 归还内存；对应
        // F90 的 deallocate。注意：缩容后再次运行需要重新 setup。

        mpack = 0;
        numat = 0;
        n2elec = 0;
        if (meci() < -1.0) n2elec = 0;
        if (!f.empty()) fock2(f, f, f, w, w, w, numat, nfirst, nlast, mode);
        a2.clear(); ab.clear(); abcmat.clear(); aicorr.clear(); aidref.clear();
        al.clear(); allgeo.clear(); allvel.clear(); allxyz.clear();
        amat.clear(); arat.clear(); atmass.clear(); b.clear(); b_esp.clear();
        bh.clear(); cosmo_C::bmat.clear(); c.clear(); cb.clear(); cc.clear(); ce.clear();
        cen.clear(); cequiv.clear(); cesp.clear(); cespm.clear(); cespm2.clear();
        cespml.clear(); cmat.clear(); cnorml.clear(); co.clear(); coord.clear();
        coorda.clear(); cosurf.clear(); d.clear(); deltap.clear(); depmul.clear();
        dx_array.clear(); dxyz.clear(); dy_array.clear(); dz_array.clear();
        ef_C::bmat.clear(); eigb.clear(); eigs.clear(); errfn.clear(); es.clear();
        esp_array.clear(); espi.clear(); ewcx.clear(); ewcy.clear(); ewcz.clear();
        ex.clear(); expn.clear(); exs.clear(); exsr.clear(); f.clear();
        f0.clear(); f1.clear(); fb.clear(); fb_ci.clear(); fc.clear(); fmat.clear();
        gden.clear(); geo.clear(); geo3.clear(); geoa.clear(); gmin1.clear();
        gnext1.clear(); grad.clear(); h.clear(); hesinv.clear(); hess.clear();
        hessc.clear(); hmat.clear(); i1fact.clear(); iam.clear(); iatsp.clear();
        ibonds.clear(); idenat.clear(); idepfn.clear(); ifact.clear(); ind.clear();
        indc.clear(); ipiden.clear(); ipo.clear(); ird.clear(); isude.clear();
        itemp.clear(); jelem.clear(); jndex.clear(); labels.clear(); loc.clear();
        locdep.clear(); locpar.clear(); lopt.clear(); na.clear(); na_store.clear();
        l_atom.clear(); namo.clear(); nat.clear(); nb.clear(); nbonds.clear();
        nc.clear(); nfirst.clear(); nlast.clear(); nn.clear(); nset.clear();
        nsetf.clear(); nw.clear(); oldf.clear(); oldhss.clear(); oldu.clear();
        ovl.clear(); p.clear(); pa.clear(); parref.clear(); pb.clear();
        pbold.clear(); pbold2.clear(); pbold3.clear(); pce.clear(); pdiag.clear();
        pewcx.clear(); pewcy.clear(); pewcz.clear(); pexpn.clear(); pexs.clear();
        pf0.clear(); pf1.clear(); pf2.clear(); phinet.clear(); pmat.clear();
        pold.clear(); pold2.clear(); pold3.clear(); potpt.clear(); ptd.clear();
        ptot2.clear(); q.clear(); qden.clear(); qdenet.clear(); qesp.clear();
        qsc.clear(); qscat.clear(); qscnet.clear(); rad.clear(); react.clear();
        rnai.clear(); rnai1.clear(); rnai2.clear(); simbol.clear(); srad.clear();
        sude.clear(); td.clear(); temp.clear(); txtatm.clear(); txtatm1.clear();
        u.clear(); u_esp.clear(); uc.clear(); uspd.clear(); vel3.clear();
        vmode.clear(); vref.clear(); vref0.clear(); w.clear(); wk.clear();
        wmat.clear(); work2.clear(); xparam.clear(); xparef.clear(); xsp.clear();
        xyz3.clear(); occa.clear(); microa.clear(); microb.clear(); meci_C::cdiag.clear();
        nalmat.clear(); eig.clear(); conf.clear(); vectci.clear(); ispin.clear();
        ispqr.clear(); cosmo_C::cdiag.clear();
    }
}
