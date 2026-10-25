// fock2.cpp — C++ translation of MOPAC 2016 "fock2.F90".
// Forms the two-electron two-center repulsion part of the Fock matrix.
// 1-based packing convention: f/ptot/p/w/wj/wk all 1-based (index 0 padding).

//
// ============================================================================
// 【中文导读】fock2() —— 用密度矩阵 P 构建 Fock 矩阵的双电子部分
// ----------------------------------------------------------------------------
// 对应 fock2.F90。iter() 每轮 SCF 都会调用：F = H + 双电子贡献(P)。
// 文件内结构：
//   - fockdorbs(): d 轨道（或大轨道块）的双中心积分贡献 (ij|kl)；
//   - fock1dorbs(): 一中心 d 轨道块；
//   - addfck(): COSMO 修正（已接线 cosmo.cpp 真实实现）；
//   - fock2(): 主入口——按原子对 (i,j) 遍历，用 W 的当前游标 kk 取积分，
//     把库仑/交换项累加进 F。分支按轻/重原子: hh（jab+kab，100 项块）、
//     hl、ll，以及周期性体系（id!=0）的四重循环。
// 数组约定：F/P/W 均为 1-based 打包下三角，索引 j*(j-1)/2+i；w 游标 kk 的
// 推进顺序必须与 Fortran 完全一致，否则 W 会错位。
// 线程模型：默认串行以保逐位一致；只有导数路径（dcart）和 OpenMP/确定性
// 并行（MOPAC_DSPAR）才使用线程私有 Fock 累加器，见文件后半说明。
// ============================================================================

#include "fock2.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>
#include "molkst_C.h"
#include "cosmo.h"
#include "cosmo_C.h"
#include "jab.h"
#include "kab.h"


namespace {

// thread_local: the derivative path (dcart) calls fock2 concurrently from
// multiple OpenMP threads; each thread keeps its own lookup tables, one-time
// init flag and per-atom density block buffer.
thread_local std::vector<int> ifact, i1fact;
thread_local int icalcn = 0;
thread_local int ione = 1;
thread_local bool lid = true;
thread_local int jindex[257] = {};
// Reused flat buffer for per-atom density blocks: row i occupies [i*82, i*82+81].
// Allocated once, grown only if a larger numat appears (MOPAC allows numat up to
// 46000, so a fixed static array would be a latent overflow).
thread_local std::vector<double> ptot2buf;
thread_local int ptot2_rows = 0;

// fockdorbs: general d-orbital (or large) two-center block.
// 中文：d 轨道（或大块）双中心贡献。ia>ja 与 ia<=ja 两种存储/遍历顺序对应
// F90 的两个分支，aa/bb 是对角块 1、非对角块 2 的权重。

void fockdorbs(int ia, int ib, int ja, int jb, std::vector<double>& f,
               const std::vector<double>& p, const std::vector<double>& ptot,
               const std::vector<double>& w, int& kr, const std::vector<int>& fact) {
    if (ia > ja) {
        for (int i = ia; i <= ib; ++i) {
            int ka = fact[i];
            double aa = 2.0;
            for (int j = ia; j <= i; ++j) {
                if (i == j) aa = 1.0;
                int kb = fact[j];
                int ij = ka + j;
                for (int k = ja; k <= jb; ++k) {
                    int kc = fact[k];
                    int ik = ka + k;
                    int jk = kb + k;
                    double bb = 2.0;
                    for (int l = ja; l <= k; ++l) {
                        if (k == l) bb = 1.0;
                        int il = ka + l;
                        int jl = kb + l;
                        int kl = kc + l;
                        kr = kr + 1;
                        double a = w[kr];
                        f[ij] += bb * a * ptot[kl];
                        f[kl] += aa * a * ptot[ij];
                        a = a * aa * bb * 0.25;
                        f[ik] -= a * p[jl];
                        f[il] -= a * p[jk];
                        f[jk] -= a * p[il];
                        f[jl] -= a * p[ik];
                    }
                }
            }
        }
    } else {
        int kref = kr;
        int nn = jb - ja + 1;
        nn = (nn * (nn + 1)) / 2;
        int n1 = 0;
        for (int i = ja; i <= jb; ++i) {
            int ka = fact[i];
            double aa = 2.0;
            for (int j = ja; j <= i; ++j) {
                n1 = n1 + 1;
                if (i == j) aa = 1.0;
                int kb = fact[j];
                int ij = ka + j;
                int n2 = 0;
                for (int k = ia; k <= ib; ++k) {
                    int kc = fact[k];
                    int ik = ka + k;
                    int jk = kb + k;
                    double bb = 2.0;
                    for (int l = ia; l <= k; ++l) {
                        n2 = n2 + 1;
                        if (k == l) bb = 1.0;
                        int il = ka + l;
                        int jl = kb + l;
                        int kl = kc + l;
                        kr = kr + 1;
                        double a = w[kref + (n2 - 1) * nn + n1];
                        f[ij] += bb * a * ptot[kl];
                        f[kl] += aa * a * ptot[ij];
                        a = a * aa * bb * 0.25;
                        f[ik] -= a * p[jl];
                        f[il] -= a * p[jk];
                        f[jk] -= a * p[il];
                        f[jl] -= a * p[ik];
                    }
                }
            }
        }
    }
}

// 2026-10 修复：删除遮蔽 cosmo.cpp 真实实现的空 addfck；调用点改为
// addfck(density, fock)（C++ 实现的参数顺序为 (pin=density, fin=Fock)）。

} // namespace

void fock1dorbs(std::vector<double>& f, const std::vector<double>& ptot,
                const std::vector<double>& pa, int mpack,
                const std::vector<double>& w,
                int& kr, int ia, int ib, int ilim) {
    // F90: w is a (ilim, ilim) array in column-major order; the caller passes
    // the base (or a kk+1 slice) of the packed w storage.  Index (ijw, klw)
    // maps to w[ijw + (klw-1)*ilim] with 1-based addressing.
    for (int i = ia; i <= ib; ++i) {
        int iw = i - ia + 1;
        for (int j = ia; j <= i; ++j) {
            int jw = j - ia + 1;
            int ij = (i * (i - 1)) / 2 + j;
            int ijw = (iw * (iw - 1)) / 2 + jw;
            double sum = 0.0;
            for (int k = ia; k <= ib; ++k) {
                int kw = k - ia + 1;
                for (int l = ia; l <= ib; ++l) {
                    int lw = l - ia + 1;
                    int ip = std::max(k, l), jp = std::min(k, l);
                    int ijp = (ip * (ip - 1)) / 2 + jp;
                    int im = std::max(kw, lw), jm = std::min(kw, lw);
                    int klw = (im * (im - 1)) / 2 + jm;
                    im = std::max(kw, jw); jm = std::min(kw, jw);
                    int ikw = (im * (im - 1)) / 2 + jm;
                    im = std::max(lw, iw); jm = std::min(lw, iw);
                    int jlw = (im * (im - 1)) / 2 + jm;
                    sum = sum + ptot[ijp] * w[ijw + (klw - 1) * ilim] -
                           pa[ijp] * w[ikw + (jlw - 1) * ilim];
                }
            }
            f[ij] += sum;
        }
    }
    kr += ilim * ilim;
}

// fock2 主入口：f 是累加器，ptot 是总密度（UHF 时为 alpha+beta）；mpack 为
// 打包长度。函数先把 ptot 复制成按原子分块的 ptot2 缓存，然后遍历原子对。

void fock2(std::vector<double>& f, const std::vector<double>& ptot,
           std::vector<double>& p, const std::vector<double>& w,
           const std::vector<double>& wj, const std::vector<double>& wk,
           int numat, const std::vector<int>& nfirst, const std::vector<int>& nlast,
           int mode) {
    using namespace molkst_C;
    bool deriv = numat < 0;
    numat = std::abs(numat);
    if (numat == 0) return;
    auto t0 = std::chrono::steady_clock::now();
    static thread_local long long n_hh = 0, n_hl = 0, n_ll = 0, n_fd = 0, n_pr = 0;
    static thread_local long long t_hh = 0, t_hl = 0, t_ll = 0, t_fd = 0, t_pr = 0;
    if (icalcn != molkst_C::numcal) {
        ifact.assign(molkst_C::norbs + 4, 0);
        i1fact.assign(molkst_C::norbs + 4, 0);
        for (int i = 1; i <= molkst_C::norbs; ++i) {
            ifact[i] = (i * (i - 1)) / 2;
            i1fact[i] = ifact[i] + i;
        }
        int m = 0;
        for (int i = 1; i <= 4; ++i)
            for (int j = 1; j <= 4; ++j) {
                int ij = std::min(i, j);
                int ji = i + j - ij;
                for (int k = 1; k <= 4; ++k) {
                    int ik = std::min(i, k);
                    for (int l = 1; l <= 4; ++l) {
                        m = m + 1;
                        int kl = std::min(k, l);
                        int lk = k + l - kl;
                        int jl = std::min(j, l);
                        jindex[m] = (ifact[ji] + ij) * 10 + ifact[lk] + kl - 10;
                    }
                }
            }
        lid = molkst_C::id == 0;
        ione = (molkst_C::id != 0) ? 0 : 1;
        icalcn = molkst_C::numcal;
    }
    // ptot2 rows are fully overwritten below (m = 1..block^2 <= 81), so no
    // clearing is needed -- avoids one heap allocation per SCF iteration.
    if (numat > ptot2_rows) {
        ptot2buf.resize((size_t)(numat + 1) * 82);
        ptot2_rows = numat;
    }
    int l = 0;
    for (int i = 1; i <= numat; ++i) {
        int ia = nfirst[i];
        int ib = nlast[i];
        int m = 0;
        for (int j = ia; j <= ib; ++j)
            for (int k = ia; k <= ib; ++k) {
                m = m + 1;
                int jk = std::min(j, k);
                int kj = k + j - jk;
                jk = jk + (kj * (kj - 1)) / 2;
                ptot2buf[(size_t)i * 82 + m] = ptot[jk];
            }
    }
    int kk = 0;
    // ---- atom-pair body shared by the serial and OpenMP paths ----

    // 以下原子对循环体在串行、OpenMP、确定性并行（MOPAC_DSPAR 分桶）三条
    // 路径间共享；kk（W 游标）随每个原子对按固定顺序推进。

    // kk (w cursor) is advanced inside; f_out is the Fock accumulator (f in
    // serial, a per-thread private copy under OpenMP).
    auto work_pair = [&](int ii, int jj, int& kk, std::vector<double>& f_out) {
        int ia = nfirst[ii];
        int ib = nlast[ii];
        int ja = nfirst[jj];
        int jb = nlast[jj];
        if (lid) {
            if (ib - ia >= 6 || jb - ja >= 6) {
                auto b0 = std::chrono::steady_clock::now();
                fockdorbs(ia, ib, ja, jb, f_out, p, ptot, w, kk, ifact);
                t_fd += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - b0).count();
                ++n_fd;
            } else if (ib - ia >= 3 && jb - ja >= 3) {
                auto b0 = std::chrono::steady_clock::now();
                // heavy-atom - heavy-atom: jab + kab, w block 100

                // 中文：重原子-重原子对，库仑项用 jab、交换项用 kab，W 每对推进 100。

                double pja[17] = {}, pjb[17] = {}, pk[17] = {};
                for (int mm = 1; mm <= 16; ++mm) {
                    pja[mm] = ptot2buf[(size_t)ii * 82 + mm];
                    pjb[mm] = ptot2buf[(size_t)jj * 82 + mm];
                }
                jab(ia, ja, pja, pjb, &w[kk], f_out.data());
                int ll = 1;
                for (int i = ia; i <= ib; ++i) {
                    int i1 = ifact[i] + ja;
                    for (int j = ll; j <= ll + 3; ++j) {
                        pk[j] = p[i1];
                        i1 = i1 + 1;
                    }
                    ll = 4 + ll;
                }
                kab(ia, ja, pk, &w[kk], f_out.data());
                kk = kk + 100;
                t_hh += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - b0).count();
                ++n_hh;
            } else if (ib - ia >= 3 && ja == jb) {
                auto b0 = std::chrono::steady_clock::now();
                // heavy-atom ii - light-atom jj  (F90 fock2.F90:190-233)
                double sumdia = 0.0, sumoff = 0.0;
                int ll = i1fact[ja];          // F90:199 ll=i1fact(ja) (light jj diagonal)
                int k = 0;
                for (int i = 0; i <= 3; ++i) {
                    int j1 = ifact[ia + i] + ia - 1;  // F90:202 j1=ifact(ia+i)+ia-1 (heavy ii)
                    if (i > 0) {
                        for (int j = 1; j <= i; ++j) {
                            f_out[j + j1] += ptot[ll] * w[j + kk + k];
                            sumoff += ptot[j + j1] * w[j + kk + k];
                        }
                        k = i + k;
                        j1 = i + j1;
                    }
                    j1 = j1 + 1;
                    k = k + 1;
                    f_out[j1] += ptot[ll] * w[kk + k];
                    sumdia += ptot[j1] * w[kk + k];
                }
                f_out[ll] = f_out[ll] + sumoff * 2.0 + sumdia;
                k = 0;
                for (int i = ia; i <= ib; ++i) {
                    int i1 = ifact[i] + ja;
                    double sum = 0.0;
                    for (int j = 1; j <= ib - ia + 1; ++j)
                        sum += p[ifact[j - 1 + ia] + ja] * w[kk + jindex[j + k]];
                    k = ib - ia + 1 + k;
                    f_out[i1] -= sum;
                }
                kk = kk + 10;
                t_hl += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - b0).count();
                ++n_hl;
            } else if (jb - ja >= 3 && ia == ib) {
                auto b0 = std::chrono::steady_clock::now();
                double sumdia = 0.0, sumoff = 0.0;
                int ll = i1fact[ia];          // F90:243 ll=i1fact(ia) (light ii diagonal)
                int k = 0;
                for (int i = 0; i <= 3; ++i) {
                    int j1 = ifact[ja + i] + ja - 1;  // F90:246 j1=ifact(ja+i)+ja-1 (heavy jj)
                    if (i > 0) {
                        for (int j = 1; j <= i; ++j) {
                            f_out[j + j1] += ptot[ll] * w[j + kk + k];
                            sumoff += ptot[j + j1] * w[j + kk + k];
                        }
                        k = i + k;
                        j1 = i + j1;
                    }
                    j1 = j1 + 1;
                    k = k + 1;
                    f_out[j1] += ptot[ll] * w[kk + k];
                    sumdia += ptot[j1] * w[kk + k];
                }
                f_out[ll] = f_out[ll] + sumoff * 2.0 + sumdia;
                k = ifact[ia] + ja;
                int j = 0;
                for (int i = k; i <= k + 3; ++i) {
                    double sum = 0.0;
                    for (int ll2 = 1; ll2 <= 4; ++ll2)
                        sum += p[ll2 - 1 + k] * w[kk + jindex[ll2 + j]];
                    j = 4 + j;
                    f_out[i] -= sum;
                }
                kk = kk + 10;
                t_hl += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - b0).count();
                ++n_hl;
            } else {
                auto b0 = std::chrono::steady_clock::now();
                // light-atom - light-atom
                int i1 = i1fact[ia];
                int j1 = i1fact[ja];
                int ij = i1 + ja - ia;
                f_out[i1] += ptot[j1] * w[kk + 1];
                f_out[j1] += ptot[i1] * w[kk + 1];
                f_out[ij] -= p[ij] * w[kk + 1];
                kk = kk + 1;
                t_ll += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - b0).count();
                ++n_ll;
            }
        } else {
            auto b0 = std::chrono::steady_clock::now();
            // periodic / solid-state general four-fold loop
            for (int i = ia; i <= ib; ++i) {
                int ka = ifact[i];
                for (int j = ia; j <= i; ++j) {
                    int kb = ifact[j];
                    int ij = ka + j;
                    double aa = 2.0;
                    if (i == j) aa = 1.0;
                    for (int k = ja; k <= jb; ++k) {
                        int kc = ifact[k];
                        int ik = (i >= k) ? ka + k : 0;
                        int jk = (j >= k) ? kb + k : 0;
                        for (int ll = ja; ll <= k; ++ll) {
                            int il = (i >= ll) ? ka + ll : 0;
                            int jl = (j >= ll) ? kb + ll : 0;
                            int kl = kc + ll;
                            double bb = 2.0;
                            if (k == ll) bb = 1.0;
                            kk = kk + 1;
                            double aj = wj[kk];
                            double ak = wk[kk];
                            if (kl > ij) continue;
                            if (i == k && aa + bb < 2.1) {
                                f_out[ij] += aj * ptot[kl];
                            } else {
                                f_out[ij] += bb * aj * ptot[kl];
                                f_out[kl] += aa * aj * ptot[ij];
                                double a = ak * aa * bb * 0.25;
                                if (jl > 0) f_out[ik] -= a * p[jl];
                                if (jk > 0) f_out[il] -= a * p[jk];
                                if (jk > 0) f_out[jk] -= a * p[il];
                                if (jl > 0) f_out[jl] -= a * p[ik];
                            }
                        }
                    }
                }
            }
            t_pr += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - b0).count();
            ++n_pr;
        }
    };
#ifdef _OPENMP

    // 并行策略：普通分子/非导数路径可选 OpenMP；论文用的确定性路径
    // （MOPAC_DSPAR）把原子对按固定顺序分到 NB 个桶、各桶私有 Fock，
    // 最后按桶序归约，保证与线程数无关。默认（无环境变量）走串行。

    // OpenMP parallel atom-pair loop (molecular, non-derivative, non-periodic,
    // mode != 2 and large enough). Each thread accumulates into a private Fock
#endif
#ifdef _OPENMP
    // Deterministic parallel atom-pair loop (MOPAC_DSPAR=1, paper path C2):
    // per-pair w-cursor (kk) offsets are computed in the serial traversal order
    // from closed-form counts (fd/hh/hl/ll/pr); pairs are dispatched to a FIXED
    // number of buckets (NB=8, independent of the OpenMP thread count); each
    // bucket accumulates a private Fock in the serial sub-order; buckets are
    // merged in fixed order 0..NB-1.  Result is therefore independent of thread
    // count / scheduling — deterministic parallel SCF.  The default serial path
    // (DSPAR unset) stays bit-identical to the official exe.
    if (!deriv && mode == 0 && std::getenv("MOPAC_DSPAR")) {
        int NB = 8;
        if (std::getenv("MOPAC_DSPAR_NB"))
            NB = std::max(1, std::min(64, std::atoi(std::getenv("MOPAC_DSPAR_NB"))));
        auto pair_kind = [&](int ii, int jj, long long& inc) -> int {
            int ia = nfirst[ii], ib = nlast[ii], ja = nfirst[jj], jb = nlast[jj];
            long long a = (long long)(ib - ia + 1) * (ib - ia + 2) / 2;
            long long b = (long long)(jb - ja + 1) * (jb - ja + 2) / 2;
            if (lid) {
                if (ib - ia >= 6 || jb - ja >= 6) { inc = a * b; return 0; } // fd
                if (ib - ia >= 3 && jb - ja >= 3) { inc = 100; return 1; }  // hh
                if (ib - ia >= 3 && ja == jb)     { inc = 10;  return 2; }  // hl
                if (jb - ja >= 3 && ia == ib)     { inc = 10;  return 2; }  // hl
                inc = 1; return 3;                                          // ll
            }
            inc = a * b; return 4;                                          // pr
        };
        // 1) scan: per-pair kk offset + kind in serial order
        std::vector<long long> kkoff, incv;
        std::vector<int> kkind, pii, pjj;
        long long kkcur = 0;
        for (int ii = 1; ii <= numat; ++ii) {
            int iminus = deriv ? ii - 1 : ii - ione;
            for (int jj = 1; jj <= iminus; ++jj) {
                long long inc = 0;
                int kd = pair_kind(ii, jj, inc);
                kkoff.push_back(kkcur); incv.push_back(inc);
                kkind.push_back(kd); pii.push_back(ii); pjj.push_back(jj);
                kkcur += inc;
            }
        }
        const size_t np = kkoff.size();
        // 2) fixed-NB buckets, round-robin dispatch by pair index (load
        //    balanced).  Bucket-internal order remains the serial sub-order,
        //    so each bucket's partial Fock is unchanged by the dispatch.
        std::vector<int> bsz(NB, 0);
        for (size_t p = 0; p < np; ++p) bsz[p % NB]++;
        std::vector<int> bstart(NB + 1, 0);
        for (int b = 0; b < NB; ++b) bstart[b + 1] = bstart[b] + bsz[b];
        std::vector<std::vector<int>> bidx(NB);
        for (size_t p = 0; p < np; ++p) bidx[p % NB].push_back((int)p);
        // per-bucket private Fock; bucket 0 carries the input f (so NB=1
        // reproduces the serial accumulation bit-identically)
        size_t fsz = f.size();
        std::vector<std::vector<double>> fbk(NB, std::vector<double>(fsz, 0.0));
        for (size_t i = 0; i < fsz; ++i) fbk[0][i] = f[i];
        // 3) process buckets in parallel.  Bucket contents depend only on the
        //    pair set inside the bucket (serial sub-order) — never on how the
        //    OpenMP runtime distributes buckets across threads — and the merge
        //    order 0..NB-1 is fixed, so the final Fock is bit-identical for
        //    ANY thread count / scheduling (deterministic parallel SCF).
#pragma omp parallel for schedule(static)
        for (int b = 0; b < NB; ++b) {
            auto& fb = fbk[b];
            for (int idx : bidx[b]) {
                int kkt = (int)kkoff[idx];
                work_pair(pii[idx], pjj[idx], kkt, fb);
            }
        }
        // 4) fixed-order merge 0..NB-1
        for (size_t i = 0; i < fsz; ++i) {
            double s = 0.0;
            for (int b = 0; b < NB; ++b) s += fbk[b][i];
            f[i] = s;
        }
        if (cosmo_C::useps) addfck(const_cast<double*>(ptot.data()), f.data());
        return;
    }
#endif
    for (int ii = 1; ii <= numat; ++ii) {
        int ia = nfirst[ii];
        int ib = nlast[ii];
        int iminus = deriv ? ii - 1 : ii - ione;
        for (int jj = 1; jj <= iminus; ++jj)
            work_pair(ii, jj, kk, f);
        if (mode == 2) {
            int ilim = ((ib - ia + 1) * (ib - ia + 2)) / 2;
            static thread_local std::vector<double> wsub;  // per-thread reuse
            if (wsub.size() < static_cast<size_t>(ilim) * ilim + 1)
                wsub.assign(static_cast<size_t>(ilim) * ilim + 1, 0.0);
            // F90 passes w(kk+1) as the (ilim,ilim) block in column-major order.
            for (int b = 1; b <= ilim; ++b)
                for (int a = 1; a <= ilim; ++a)
                    wsub[a + (b - 1) * ilim] = w[kk + (a - 1) + (b - 1) * ilim + 1];
            fock1dorbs(f, ptot, p, molkst_C::mpack, wsub, kk, ia, ib, ilim);
        }
    }
    if (cosmo_C::useps) addfck(const_cast<double*>(ptot.data()), f.data());
}
