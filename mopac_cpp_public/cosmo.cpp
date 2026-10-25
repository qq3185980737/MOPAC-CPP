// ============================================================================
// 【中文导读】cosmo.cpp —— COSMO 溶剂化模型（cavity / B 矩阵 / 介电求解）
// ----------------------------------------------------------------------------
// 对应 cosmo.F90。主要流程：
//   1) cosini()  初始化半径 srad、原子对列表 nipa/lipa、方向向量 dirvec/dirsm、
//                n0、disex2 等；
//   2) coscav()  构造 SAS 腔体：逐原子撒点、按 Segments 归属、surclo() 闭合凹面
//                环段、累加 area/cosvol，并用 mfinel() 生成 finel 列表、组装并
//                做 Cholesky 分解得到 amat；
//   3) mkbmat()  B 矩阵（网格/环段到原子的 1/r 多极作用）；
//   4) addnuc()/addfck() 计算核/电子屏蔽电荷 qscnet 与介电能量 ediel，并把
//                B*qscel 反馈到 Fock 矩阵；
//   5) coscl1()/coscl2() 打包下三角 Cholesky 分解与回代求解。
// 关键约定：rsc、dirvec、dirsm、finel 等按 Fortran 列主序打包；cosurf/iatsp/
// nsetf 等为 1-based。COSMO 的数值闭环已在 2026-10 与官方逐项对照（见
// docs/TRANSLATION_GAPS.md §19-§26）。
// ============================================================================


// cosmo.cpp 鈥?C++ translation of MOPAC 2016 "cosmo.F90".
// 2026-10: coscav/surclo/mkbmat/cosini 已完整移植（曾为占位）。

#include "cosmo.h"

#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>

#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "mopend.h"
#include "parameters_C.h"
#include "reada.h"

void mxm(const double* a, int lda, const double* b, int n, double* c, int ldc);

using namespace common_arrays_C;
using namespace cosmo_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace parameters_C;

// Packed lower-triangle Cholesky factorization: A = L*L^T, stored in A.
// id(i) = i(i-1)/2, Fortran 1-based. info<0 signals non-positive pivot.
void coscl1(double* a, std::vector<int>& id, int n, int& info) {
    id.resize(n + 1);
    int indi = 0;
    for (int i = 1; i <= n; ++i) { id[i] = indi; indi += i; }
    info = 0;
    for (int k = 1; k <= n; ++k) {
        int indk = id[k];
        int kk = k + indk;
        for (int i = k; i <= n; ++i) {
            indi = id[i];
            double summe = 0.0;
            for (int j = 1; j <= k - 1; ++j)
                summe += a[j + indi] * a[j + indk];
            summe = a[k + indi] - summe;
            if (i == k) {
                if (summe < 0.0) { info = -1; summe = a[kk]; }
                a[kk] = 1.0 / std::sqrt(summe);
            } else {
                a[k + indi] = summe * a[kk];
            }
        }
    }
}

// Solve C x = y given the Cholesky factors in A (packed).
void coscl2(const double* a, const std::vector<int>& id,
            double* x, const double* y, int n) {
    for (int k = 1; k <= n; ++k) {
        double summe = y[k];
        int kk = id[k];
        for (int i = k - 1; i >= 1; --i)
            summe -= a[i + kk] * x[i];
        x[k] = summe * a[k + kk];
    }
    for (int k = n; k >= 1; --k) {
        double summe = x[k];
        for (int i = k + 1; i <= n; ++i)
            summe -= a[k + id[i]] * x[i];
        x[k] = summe * a[k + id[k]];
    }
}

namespace { void coscl2_dummy() {} }

// CI-state dielectric correction: build CI density difference, back-transform
// to AO basis, solve screening charges, add dielectric energy to eig.
void dmecip(double* coeffs, double* deltap, double* delta, double* eig,
            double* vectci, const double* conf) {
    using namespace meci_C;
    double fact = -ev * a0 * fnsq / 2.0;
    for (int ist = 1; ist <= lab; ++ist) {
        int ii = lab * (ist - 1);
        for (int j = 1; j <= lab; ++j) vectci[j] = conf[ii + j];
        for (int i = 1; i <= nmos; ++i) {
            deltap[(i-1)*nmos + i-1] = -occa[i] * 2.0;
            for (int j = 1; j <= i - 1; ++j)
                deltap[(i-1)*nmos + j-1] = 0.0;
        }
        for (int id = 1; id <= lab; ++id)
            for (int jd = 1; jd <= id; ++jd) {
                if (nalmat[id] != nalmat[jd]) continue;
                int ix = 0, iy = 0;
                for (int j = 1; j <= nmos; ++j) {
                    ix += std::abs(microa[j][id]-microa[j][jd]);
                    iy += std::abs(microb[j][id]-microb[j][jd]);
                }
                if (ix + iy > 2) continue;
                if (ix == 2) {
                    int i; for (i = 1; i <= nmos; ++i)
                        if (microa[i][id] != microa[i][jd]) break;
                    int ij = microb[i][id];
                    int j;
                    for (j = i+1; j <= nmos; ++j) {
                        if (microa[j][id] != microa[j][jd]) break;
                        ij += microa[j][id] + microb[j][id];
                    }
                    double sum = 0;
                    for (int k = 1; k <= nstate; ++k)
                        sum += vectci[id+(k-1)*lab]*vectci[jd+(k-1)*lab];
                    deltap[(j-1)*nmos + i-1] += sum * (1 - 2*(ij%2)) / nstate;
                } else if (iy == 2) {
                    int i; for (i = 1; i <= nmos; ++i)
                        if (microb[i][id] != microb[i][jd]) break;
                    int ij = 0;
                    int j; for (j = i+1; j <= nmos; ++j) {
                        if (microb[j][id] != microb[j][jd]) break;
                        ij += microa[j][id] + microb[j][id];
                    }
                    ij += microa[j][id];
                    double sum = 0;
                    for (int k = 1; k <= nstate; ++k)
                        sum += vectci[id+(k-1)*lab]*vectci[jd+(k-1)*lab];
                    deltap[(j-1)*nmos + i-1] += sum * (1 - 2*(ij%2)) / nstate;
                } else {
                    double sum = 0;
                    for (int k = 1; k <= nstate; ++k)
                        sum += vectci[id+(k-1)*lab]*vectci[id+(k-1)*lab];
                    for (int i = 1; i <= nmos; ++i)
                        deltap[(i-1)*nmos + i-1] += (microa[i][id]+microb[i][id])*sum/nstate;
                }
            }
        for (int i = 1; i <= nmos; ++i)
            for (int j = 1; j <= i - 1; ++j)
                deltap[(j-1)*nmos + i-1] = deltap[(i-1)*nmos + j-1];
        // delta = coeffs(1:norbs, nelec+1:nelec+nmos) * deltap
        mxm(coeffs + (nelec), norbs, deltap, nmos, delta, nmos);
        int iden = 0;
        for (int iat = 1; iat <= numat; ++iat)
            for (int i = nfirst[iat]; i <= nlast[iat]; ++i)
                for (int j = nfirst[iat]; j <= i; ++j) {
                    double sum = 0;
                    for (int k = 1; k <= nmos; ++k)
                        sum += delta[(i-1)*nmos + k-1] * coeffs[(j-1)*norbs + nelec + k-1];
                    ++iden;
                    qdenet[iden][1] = sum;
                }
        for (int i = 1; i <= nps; ++i) {
            double sum = 0;
            for (int j = 1; j <= lm61; ++j) sum += bmat[j][i]*qdenet[j][1];
            phinet[i][1] = sum;
        }
        // coscl2 闇€瑕佽繛缁?1..nps 鐨勪竴鍒楋紱C++ 瀛樺偍鏄?[segment][slot]锛?
        // 鍥犳鍏堝鍒?slot-1 鍒楀埌鎵佸钩缂撳啿锛岃В瀹屽啀鍐欏洖锛堝悓 addnuc锛夈€?
        {
            std::vector<double> q_flat(nps + 1, 0.0), p_flat(nps + 1, 0.0);
            std::vector<int> id_c2(nps + 1, 0);
            for (int k = 1; k <= nps; ++k) p_flat[k] = phinet[k][1];
            for (int i = 2; i <= nps; ++i) id_c2[i] = id_c2[i - 1] + i - 1;
            coscl2(amat.data(), id_c2, &q_flat[0], &p_flat[0], nps);
            for (int k = 1; k <= nps; ++k) qscnet[k][1] = q_flat[k];
        }
        double edie = 0;
        for (int i = 1; i <= nps; ++i) edie += qscnet[i][1]*phinet[i][1];
        eig[ist] += fact * edie;
    }
    for (int j = 1; j <= lab; ++j) vectci[j] = conf[j];
}

// COSMO analytic gradient contribution (electrostatic + surface closure).
void diegrd(double* dxyz) {
    double db[4][11];
    for (int i = 1; i <= 10; ++i)
        for (int ix = 1; ix <= 3; ++ix) db[ix][i] = 0.0;
    db[0][1] = 1.0;
    double fact = -ev * a0 * fpc_9;
    for (int k = 1; k <= nps; ++k) {
        int iak = iatsp[k];
        double xk1 = cosurf[1][k], xk2 = cosurf[2][k], xk3 = cosurf[3][k];
        double qsk = qscnet[k][3];
        for (int l = 1; l <= k - 1; ++l) {
            int ial = iatsp[l];
            if (ial != iak) {
                double xl1 = cosurf[1][l] - xk1;
                double xl2 = cosurf[2][l] - xk2;
                double xl3 = cosurf[3][l] - xk3;
                double d2 = xl1*xl1 + xl2*xl2 + xl3*xl3;
                double ff = qsk * qscnet[l][3] * fact * std::pow(d2, -1.5) / fepsi;
                dxyz[(iak-1)*3+0] -= xl1*ff; dxyz[(ial-1)*3+0] += xl1*ff;
                dxyz[(iak-1)*3+1] -= xl2*ff; dxyz[(ial-1)*3+1] += xl2*ff;
                dxyz[(iak-1)*3+2] -= xl3*ff; dxyz[(ial-1)*3+2] += xl3*ff;
            }
        }
    }
    double bsurf = 0.0;
    for (int i = 1; i <= nipc; ++i) {
        int ia = isude[1][i], ib = isude[2][i];
        double deab = -0.25 * (qscat[ia]*qscat[ia]*sude[1][i]/arat[ia]
                             + qscat[ib]*qscat[ib]*sude[2][i]/arat[ib]
                             + bsurf*(sude[1][i]+sude[2][i]));
        double x1 = coord[0][ib]-coord[0][ia];
        double x2 = coord[1][ib]-coord[1][ia];
        double x3 = coord[2][ib]-coord[2][ia];
        deab /= std::sqrt(x1*x1+x2*x2+x3*x3);
        dxyz[(ia-1)*3+0] -= x1*deab; dxyz[(ib-1)*3+0] += x1*deab;
        dxyz[(ia-1)*3+1] -= x2*deab; dxyz[(ib-1)*3+1] += x2*deab;
        dxyz[(ia-1)*3+2] -= x3*deab; dxyz[(ib-1)*3+2] += x3*deab;
    }
    for (int k = 1; k <= nps; ++k) {
        int iak = iatsp[k];
        double xk1 = cosurf[1][k], xk2 = cosurf[2][k], xk3 = cosurf[3][k];
        double qsk = qscnet[k][3];
        int iden = 0;
        for (int i = 1; i <= numat; ++i) {
            int idel = nlast[i]+1 - nfirst[i];
            if (i != iak) {
                int nati = nat[i];
                double xx1 = xk1-coord[0][i], xx2 = xk2-coord[1][i], xx3 = xk3-coord[2][i];
                double d2 = xx1*xx1+xx2*xx2+xx3*xx3;
                double ddi = dd[nati]*a0;
                double qqi2 = (a0*qq[nati])*(a0*qq[nati]);
                double ff0 = -qsk*fact*std::pow(d2,-1.5);
                if (idel > 1) {
                    double rm2 = 1.0/d2, rm4 = rm2*rm2;
                    db[0][2]=ddi*3*xx1*rm2; db[0][4]=ddi*3*xx2*rm2; db[0][7]=ddi*3*xx3*rm2;
                    db[0][3]=1.0+qqi2*(15*xx1*xx1*rm2-3.0)*rm2;
                    db[0][6]=1.0+qqi2*(15*xx2*xx2*rm2-3.0)*rm2;
                    db[0][10]=1.0+qqi2*(15*xx3*xx3*rm2-3.0)*rm2;
                    db[0][5]=qqi2*15*xx1*xx2*rm4;
                    db[0][8]=qqi2*15*xx1*xx3*rm4;
                    db[0][9]=qqi2*15*xx3*xx2*rm4;
                    db[1][2]=ddi; db[2][4]=db[1][2]; db[3][7]=db[1][2];
                    db[1][3]=6*qqi2*xx1*rm2;
                    db[2][6]=6*qqi2*xx2*rm2;
                    db[3][10]=6*qqi2*xx3*rm2;
                    db[1][5]=db[2][6]; db[2][5]=db[1][3];
                    db[1][8]=db[3][10]; db[3][8]=db[1][3];
                    db[2][9]=db[3][10]; db[3][9]=db[2][6];
                }
                int nj = std::min(10, idel*(idel+1)/2);
                for (int j = 1; j <= nj; ++j) {
                    double ff = -ff0*qdenet[iden+j][3];
                    if (j==1 && idel==9)
                        for (int iii=5;iii<=9;++iii)
                            ff -= ff0*qdenet[iden+(iii*(iii+1))/2][3];
                    double dx1=(xx1*db[0][j]-db[1][j])*ff;
                    double dx2=(xx2*db[0][j]-db[2][j])*ff;
                    double dx3=(xx3*db[0][j]-db[3][j])*ff;
                    dxyz[(iak-1)*3+0]+=dx1; dxyz[(i-1)*3+0]-=dx1;
                    dxyz[(iak-1)*3+1]+=dx2; dxyz[(i-1)*3+1]-=dx2;
                    dxyz[(iak-1)*3+2]+=dx3; dxyz[(i-1)*3+2]-=dx3;
                }
            }
            iden += idel*(idel+1)/2;
        }
    }
}

void ansude(double ra, double rb, double d, double rs,
            double& ara, double& arb, double& aar, double& abr,
            double& arad, double& arbd, double& rinc) {
    double qa = ra + rs, qb = rb + rs;
    double ca = (qa * qa + d * d - qb * qb) / (2.0 * qa * d);
    double cbv = (qb * qb + d * d - qa * qa) / (2.0 * qb * d);
    double sa = std::sqrt(1.0 - ca * ca);
    double sb = std::sqrt(1.0 - cbv * cbv);
    double ta = pi * sa, tb = pi * sb;
    double fza = (1.0 - std::cos(ta)) / 2.0;
    double fzb = (1.0 - std::cos(tb)) / 2.0;
    if (sa < 0 || sb < 0) fza = 1.0;
    if (sa < 0 || sb < 0) fzb = 1.0;
    double xa = fzb * rs * (ca + cbv);
    double xb = fza * rs * (ca + cbv);
    double ya = ra * sa - fzb * rb * sb;
    double yb = rb * sb - fza * ra * sa;
    double za = std::sqrt(xa * xa + ya * ya);
    double zb = std::sqrt(xb * xb + yb * yb);
    rinc = 0.5 * (za + zb) / std::sqrt(rs * rs * (ca + cbv) * (ca + cbv) +
                                        (ra * sa - rb * sb) * (ra * sa - rb * sb));
    ara = pi * ra * (2.0 * (1.0 + ca) * ra + sa * za);
    arb = pi * rb * (2.0 * (1.0 + cbv) * rb + sb * zb);
    aar = pi * ra * (sa * za);
    abr = pi * rb * (sb * zb);
    double cad = (qb * qb + d * d - qa * qa) / (2.0 * qa * d * d);
    double cbd = (qa * qa + d * d - qb * qb) / (2.0 * qb * d * d);
    double sad = -ca * cad / sa, sbd = -cbv * cbd / sb;
    double tad = pi * sad, tbd = pi * sbd;
    double fzad = std::sin(ta) * 0.5, fzbd = std::sin(tb) * 0.5;
    if (sa < 0 || sb < 0) fzad = 0.0;
    if (sa < 0 || sb < 0) fzbd = 0.0;
    double xad = rs * ((ca + cbv) * fzbd * tbd + fzb * (cad + cbd));
    double xbd = rs * ((ca + cbv) * fzad * tad + fza * (cad + cbd));
    double yad = ra * sad - fzbd * tbd * rb * sb - fzb * rb * sbd;
    double ybd = rb * sbd - fzad * tad * ra * sa - fza * ra * sad;
    double zad = (xa * xad + ya * yad) / za;
    double zbd = (xb * xbd + yb * ybd) / zb;
    arad = pi * ra * (sad * za + sa * zad + 2.0 * ra * cad);
    arbd = pi * rb * (sbd * zb + sb * zbd + 2.0 * rb * cbd);
}

void ciint(const double* c34, double* pq34) {
    int i0 = 0;
    for (int i = 1; i <= lm61; ++i) pq34[i] = 0.0;
    if (nps < 0) return;
    for (int i = 1; i <= lm61; ++i) {
        for (int j = 1; j < i; ++j) {
            i0++;
            pq34[j] += cmat[i0] * c34[i];
            pq34[i] += cmat[i0] * c34[j];
        }
        i0++;
        pq34[i] += cmat[i0] * c34[i];
    }
}

void addfck(double* pin, double* fin) {
    double fcon = a0 * ev;
    for (int i = 1; i <= numat; ++i) qscat[i] = 0.0;
    for (int i = 1; i <= lm61; ++i) {
        qdenet[i][2] = gden[i] * pin[ipiden[i]];
        qdenet[i][3] = qdenet[i][2] + qdenet[i][1];
    }
    for (int i = 1; i <= nps; ++i) {
        double phi = 0.0;
        for (int j = 1; j <= lm61; ++j) phi += bmat[j][i] * qdenet[j][2];
        phinet[i][2] = phi;
        phinet[i][3] = phinet[i][1] + phi;
    }
    // coscl2 闇€瑕佽繛缁?1..nps 鐨勪竴鍒楋紱C++ 瀛樺偍鏄?[segment][slot]锛岀洿鎺ュ彇
    // &qscnet[1][2] 浼氳秺鐣屽啓鍧忓唴瀛橈紙COSMO = NaN 鐨勬牴鍥狅級锛屾晠鐢ㄦ墎骞崇紦鍐层€?
    {
        std::vector<double> q_flat(nps + 1, 0.0), p_flat(nps + 1, 0.0);
        std::vector<int> id_c2(nps + 1, 0);
        for (int k = 1; k <= nps; ++k) p_flat[k] = phinet[k][2];
        for (int i = 2; i <= nps; ++i) id_c2[i] = id_c2[i - 1] + i - 1;
        coscl2(amat.data(), id_c2, &q_flat[0], &p_flat[0], nps);
        for (int k = 1; k <= nps; ++k) qscnet[k][2] = q_flat[k];
    }
    ediel = 0.0;
    double s1 = 0.0, s3 = 0.0;
    for (int i = 1; i <= nps; ++i) {
        int iat = iatsp[i];
        qscnet[i][2] = -fepsi * qscnet[i][2];
        double qsc3 = qscnet[i][1] + qscnet[i][2];
        qscnet[i][3] = qsc3;
        ediel += qsc3 * phinet[i][3];
        s1 += qscnet[i][1];
        s3 += qsc3;
        qscat[iat] += qsc3;
    }
    ediel *= fcon / 2.0;
    for (int i = 1; i <= lm61; ++i) {
        int im = ipiden[i];
        double fim = 0.0;
        for (int j = 1; j <= nps; ++j) fim += bmat[i][j] * qscnet[j][2];
        fin[im] -= fcon * fim;
    }

}

void addhcr() {
    double fcon = a0 * ev;
    for (int i = 1; i <= lm61; ++i) {
        int im = ipiden[i];
        double him = 0.0;
        for (int j = 1; j <= nps; ++j) him += bmat[i][j] * qscnet[j][1];
        h[im] -= fcon * him;
    }
}

// Build the screening interaction matrix B (1/r Coulomb multipoles).
// Build the list of basic grid points + ring segments for segment ips.
// finel is (4, 300, 2); rsc is (4, maxrs); dirvec (4,1082); tm (3,3).
void mfinel(int ips, int k, double* finel,
            const std::vector<int>& nar_csm, const std::vector<int>& nsetf,
            const std::vector<int>& nset, double* rsc,
            const std::vector<int>& nipsrs,
            double* dirvec, double* tm, const double* x, double r,
            int& nfl, int ioldcv, int lenabc, int maxrs) {
    auto F=[&](int ix,int j,int kk)->double& { return finel[((kk-1)*300 + (j-1))*4 + (ix-1)]; };
    auto DV=[&](int ix,int idir)->double& { return dirvec[(idir-1)*4+(ix-1)]; };
    auto TM=[&](int ix,int j)->double& { return tm[(j-1)*3+(ix-1)]; };
    auto RS=[&](int ix,int irs)->double& { return rsc[(irs-1)*4+(ix-1)]; };
    nfl = 0;
    int nari = nar_csm[ips];
    for (int l = nsetf[ips]; l <= nsetf[ips]+nari-1; ++l) {
        int idir = nset[l];
        ++nfl;
        double y1=DV(1,idir)*r, y2=DV(2,idir)*r, y3=DV(3,idir)*r;
        F(1,nfl,k)=y1*TM(1,1)+y2*TM(2,1)+y3*TM(3,1)+x[0];
        F(2,nfl,k)=y1*TM(1,2)+y2*TM(2,2)+y3*TM(3,2)+x[1];
        F(3,nfl,k)=y1*TM(1,3)+y2*TM(2,3)+y3*TM(3,3)+x[2];
        F(4,nfl,k)=DV(4,idir)*r*r;
    }
    if (ioldcv == 1) return;
    int irs0 = (ips > 1) ? nipsrs[ips-1] + 1 : 1;
    int irs1 = nipsrs[ips];
    for (int irs = irs0; irs <= irs1; ++irs) {
        ++nfl;
        for (int ix = 1; ix <= 4; ++ix) F(ix,nfl,k) = RS(ix,irs);
    }
}

// Construct/update the solvent-accessible surface (SAS) cavity tessellation.
void coscav() {
    int maxrs = 70 * numat;
    auto DD=[&](int ix,int i)->double& { return (&dirvec[0][0])[(i-1)*4+(ix-1)]; };
    auto DS=[&](int ix,int i)->double& { return (&dirsm[0][0])[(i-1)*4+(ix-1)]; };
    std::vector<double> rdat(numat + 1);
    std::vector<double> rsc(4 * (maxrs + 1), 0.0);  // 鍒椾富搴? (i-1)*4+ix
    std::vector<std::array<std::array<double,3>,3>> Tmat(numat + 1);
    std::vector<std::array<int,3>> nn(numat + 1);
    std::vector<int> isort(maxrs + 1), ipsrs(maxrs + 1);
    std::vector<int> nipsrs(lenabc + 1);
    std::vector<int> nset(1082 * numat + 1);
    std::vector<int> nipa(numat + 1);
    std::vector<double> coord_flat, finel_flat, dirvec_flat;
    std::vector<double*> rsc_ptrs;
    std::vector<char> din;
    std::vector<std::array<double,3>> dirtm(1083);
    std::vector<std::array<std::array<std::array<double,4>,2>,300>> finel;
    finel.resize(301);
    if (srad.empty()) srad.resize(numat + 1);
    if (nar_csm.empty()) nar_csm.resize(lenabc + 2);

    // ---- coscav differential dump (env COSCAV_DUMP=<prefix>) ----
    std::string coscav_dump;
    if (const char* dp = std::getenv("COSCAV_DUMP")) {
        coscav_dump = dp;
        std::ofstream f(coscav_dump + ".in");
        f << std::setprecision(17);
        auto W = [&](double v) { f << v << '\n'; };
        W(numat); W(lenabc); W((double)ioldcv); W(rsolv); W(disex2);
        W((double)n0[1]); W((double)n0[2]);
        for (int idir = 1; idir <= 1082; ++idir)
            for (int ix = 1; ix <= 4; ++ix) W(DD(ix, idir));
        int nds = n0[1] + n0[2];
        for (int id = 1; id <= nds; ++id)
            for (int ix = 1; ix <= 4; ++ix) W(DS(ix, id));
        for (int a = 1; a <= numat; ++a)
            for (int ix = 0; ix < 3; ++ix) W(coord[ix][a]);
        for (int a = 1; a <= numat; ++a) W(srad[a]);
        for (int a = 1; a <= numat; ++a) W((double)nat[a]);
    }

    // MAKE COORDINATES A BIT ASYMMETRIC IN ORDER TO AVOID
    // SYMMETRY PROBLEMS WITH CAVITY CONSTRUCTION (F90 cosmo.F90 235-238)
    for (int i = 1; i <= numat; ++i)
        for (int j = 1; j <= 3; ++j)
            coord[j-1][i] += std::cos(i*j*0.1) * 3.0e-9;

    int ilipa = 0;
    for (int i = 1; i <= numat; ++i) {
        double ri = srad[i], r = ri + rsolv, rr = r + rsolv;
        double xa1 = coord[0][i], xa2 = coord[1][i], xa3 = coord[2][i];
        for (int j = 1; j <= numat; ++j) {
            if (j == i) continue;
            double d1 = xa1-coord[0][j], d2 = xa2-coord[1][j], d3 = xa3-coord[2][j];
            double dist = d1*d1+d2*d2+d3*d3;
            if (dist < (rr+srad[j])*(rr+srad[j])) ++ilipa;
        }
    }
    // lipa is addressed 1-based (lipa[ilipa]); allocate one spare slot.
    std::vector<int> lipa(std::max(1, ilipa) + 1);
    double fdiagr = 2.1 * std::sqrt(pi);
    int inset = 1;
    ilipa = 0;
    nps = 0; area = 0; cosvol = 0;

    for (int i = 1; i <= numat; ++i) {
        nipa[i] = 0;
        double ri = srad[i], r = ri + rsolv, rr = r + rsolv, ri2 = ri*ri;
        double xa1 = coord[0][i], xa2 = coord[1][i], xa3 = coord[2][i];
        int nps0 = nps + 1;
        for (int j = 1; j <= numat; ++j) {
            if (j == i) continue;
            double d1=xa1-coord[0][j], d2=xa2-coord[1][j], d3=xa3-coord[2][j];
            double dist = d1*d1+d2*d2+d3*d3;
            if (dist < (rr+srad[j])*(rr+srad[j])) {
                ++ilipa;
                if (ilipa > maxrs) { mopend("Solvent radius too large - reduce RSOLV"); return; }
                ++nipa[i]; lipa[ilipa] = j;
            }
        }
        double dist1=1e20, dist2=1e20, dist3=1e20;
        int nn1=0, nn2=0, nn3=0;
        for (int j = 1; j <= numat; ++j) {
            if (j == i) continue;
            double d1=xa1-coord[0][j], d2=xa2-coord[1][j], d3=xa3-coord[2][j];
            double dist = d1*d1+d2*d2+d3*d3;
            if (dist+0.05 < dist3) { dist3 = dist; nn3 = j; }
            if (dist3+0.05 < dist2) { double t=dist2; dist2=dist3; dist3=t; nn3=nn2; nn2=j; }
            if (dist2+0.05 < dist1) { double t=dist1; dist1=dist2; dist2=t; nn2=nn1; nn1=j; }
        }
        auto& TM = Tmat[i];
        if (nn1 == 0) { TM[0][0]=1; TM[0][1]=0; TM[0][2]=0; }
        else {
            double d1=xa1-coord[0][nn1], d2=xa2-coord[1][nn1], d3=xa3-coord[2][nn1];
            double dist = 1.0/std::sqrt(d1*d1+d2*d2+d3*d3);
            TM[0][0]=(coord[0][nn1]-xa1)*dist;
            TM[0][1]=(coord[1][nn1]-xa2)*dist;
            TM[0][2]=(coord[2][nn1]-xa3)*dist;
        }
        bool done2 = false;
        while (!done2) {
            if (nn2 == 0) {
                if (std::abs(TM[0][1])+std::abs(TM[0][0]) > 0.1) {
                    TM[1][0]=-TM[0][1]; TM[1][1]=TM[0][0]; TM[1][2]=0;
                } else {
                    TM[1][0]=-TM[0][1]; TM[1][1]=TM[0][2]; TM[1][2]=0;
                }
                done2 = true;
            } else {
                double d1=xa1-coord[0][nn2], d2=xa2-coord[1][nn2], d3=xa3-coord[2][nn2];
                double dist = 1.0/std::sqrt(d1*d1+d2*d2+d3*d3);
                double xx1=(coord[0][nn2]-xa1)*dist, xx2=(coord[1][nn2]-xa2)*dist, xx3=(coord[2][nn2]-xa3)*dist;
                double sp = xx1*TM[0][0]+xx2*TM[0][1]+xx3*TM[0][2];
                if (sp*sp > 0.99) { nn2=nn3; nn3=0; }
                else {
                    double si = 1.0/std::sqrt(1.0-sp*sp);
                    TM[1][0]=(xx1-sp*TM[0][0])*si;
                    TM[1][1]=(xx2-sp*TM[0][1])*si;
                    TM[1][2]=(xx3-sp*TM[0][2])*si;
                    done2 = true;
                }
            }
        }
        TM[2][0]=TM[0][1]*TM[1][2]-TM[1][1]*TM[0][2];
        TM[2][1]=TM[0][2]*TM[1][0]-TM[1][2]*TM[0][0];
        TM[2][2]=TM[0][0]*TM[1][1]-TM[1][0]*TM[0][1];
        for (int j = 1; j <= 1082; ++j) {
            double xx1=DD(1,j), xx2=DD(2,j), xx3=DD(3,j);
            for (int ix = 0; ix < 3; ++ix)
                dirtm[j][ix] = xx1*TM[0][ix]+xx2*TM[1][ix]+xx3*TM[2][ix];
        }
        int narea = 0;
        din.assign(1083, 0);
        std::vector<int> iseg(1083, 0);
        for (int j = 1; j <= 1082; ++j) {
            double xx1=xa1+dirtm[j][0]*r, xx2=xa2+dirtm[j][1]*r, xx3=xa3+dirtm[j][2]*r;
            bool inside = false;
            for (int ik = ilipa-nipa[i]+1; ik <= ilipa; ++ik) {
                int k = lipa[ik];
                double e1=xx1-coord[0][k], e2=xx2-coord[1][k], e3=xx3-coord[2][k];
                double dist = std::sqrt(e1*e1+e2*e2+e3*e3)-rsolv-srad[k];
                if (dist < 0) { inside = true; break; }
            }
            if (inside) continue;
            ++narea;
            cosvol += ri2*DD(4,j)*(dirtm[j][0]*xa1+dirtm[j][1]*xa2+dirtm[j][2]*xa3+ri);
            area += ri2*DD(4,j);
            din[j] = 1;
        }
        if (narea != 0) {
            int i0 = (nat[i]==1) ? 2 : 1;
            int jmax = n0[i0];
            int base = (i0-1)*n0[1];
            for (int j = 1; j <= jmax; ++j) {
                ++nps;
                if (nps > lenabc) { mopend("NPS IS GREATER THAN LENABC-USE SMALLER NSPA"); goto done_coscav; }
                iatsp[nps] = i;
                double xx1=DS(1,base+j), xx2=DS(2,base+j), xx3=DS(3,base+j);
                for (int ix = 0; ix < 3; ++ix)
                    cosurf[ix+1][nps] = xx1*TM[0][ix]+xx2*TM[1][ix]+xx3*TM[2][ix];
            }
            int niter = 0;
            do {
                ++niter;
                for (int ips = nps0; ips <= nps; ++ips) {
                    nar_csm[ips]=0;
                    for (int ix=1;ix<=3;++ix) phinet[ips][ix]=0;
                }
                for (int j = 1; j <= 1082; ++j) {
                    if (!din[j]) continue;
                    double spm=-1, x1=dirtm[j][0], x2=dirtm[j][1], x3=dirtm[j][2];
                    int ipm = 0;
                    for (int ips=nps0; ips<=nps; ++ips) {
                        double sp = x1*cosurf[1][ips]+x2*cosurf[2][ips]+x3*cosurf[3][ips];
                        if (sp >= spm) { spm = sp*(1+1e-14); ipm = ips; }
                    }
                    iseg[j]=ipm;
                    ++nar_csm[ipm];
                    for (int ix=0;ix<3;++ix)
                        phinet[ipm][ix+1] += dirtm[j][ix]*DD(4,j);
                }
                int ips = nps0-1;
                while (true) {
                    ++ips;
                    while (nar_csm[ips]==0) {
                        niter=1; --nps;
                        if (ips > nps) break;
                        for (int jps=ips;jps<=nps;++jps) {
                            nar_csm[jps]=nar_csm[jps+1];
                            for (int ix=1;ix<=3;++ix) phinet[jps][ix]=phinet[jps+1][ix];
                        }
                    }
                    if (ips > nps) break;
                    double ds2=phinet[ips][1]*phinet[ips][1]+phinet[ips][2]*phinet[ips][2]+phinet[ips][3]*phinet[ips][3];
                    ds2 = std::max(ds2,1e-20);
                    double dist = 1.0/std::sqrt(ds2);
                    for (int ix=1;ix<=3;++ix) cosurf[ix][ips] = phinet[ips][ix]*dist;
                    if (ips >= nps) break;
                }
                if (niter >= 2) break;
            } while (true);
            for (int ips=nps0; ips<=nps; ++ips) {
                nsetf[ips]=inset; inset += nar_csm[ips]; nar_csm[ips]=0;
                cosurf[4][ips]=0;
                for (int ix=0;ix<3;++ix) cosurf[ix+1][ips]=cosurf[ix+1][ips]*ri+coord[ix][i];
            }
            for (int j=1;j<=1082;++j) {
                if (!din[j]) continue;
                int ipm=iseg[j];
                int nara=nar_csm[ipm];
                nset[nsetf[ipm]+nara]=j;
                nar_csm[ipm]=nara+1;
                cosurf[4][ipm]+=DD(4,j)*ri2;
            }
        }
        if (ilipa >= 70*numat) break;
    }
    if (ioldcv == 0) {
        // Fortran cosmo.F90 L548-553: din(1..numat) is reused as an atom
        // marker before surclo: all true, then false for atoms that carry
        // at least one segment (iatsp).  surclo skips atoms with din==true.
        for (int i = 1; i <= numat; ++i) din[i] = 1;
        for (int j = 1; j <= nps; ++j) din[iatsp[j]] = 0;
        // call surclo 鈥?closure of the concave regions of the cavity.
        // C++ surclo expects a planar [3][numat+1] coordinate buffer and a
        // bool[] din; convert from the row-major global coord and vector<char>.
        coord_flat.assign((numat + 1) * 3, 0.0);
        for (int a = 1; a <= numat; ++a)
            for (int ix = 0; ix < 3; ++ix)
                coord_flat[ix * (numat + 1) + a] = coord[ix][a];
        surclo(coord_flat.data(), nipa.data(), lipa.data(),
               reinterpret_cast<const bool*>(din.data()),
               std::max(numat, 1082), rsc.data(), isort.data(),
               ipsrs.data(), nipsrs.data(), nat.data(), srad.data(), maxrs);
    }
    cosvol /= 3;
    for (int i=1;i<=numat;++i) { arat[i]=0; rdat[i]=0; }
    // FILLING AMAT (F90 cosmo.F90 lines 554-610): diagonal via mfinel
    // self-energy terms, off-diagonal via 1/distance (or mfinel pair sums
    // inside the exclusion sphere); then Cholesky factorization (coscl1).
    finel_flat.assign(2 * 300 * 4, 0.0);
    dirvec_flat.assign((1082 + 1) * 4, 0.0);
    for (int id = 1; id <= 1082; ++id)
        for (int ix = 1; ix <= 4; ++ix) dirvec_flat[(id - 1) * 4 + (ix - 1)] = DD(ix, id);
    for (int ips=1; ips<=nps; ++ips) {
        int i=iatsp[ips];
        double ri = srad[i];
        double xi[3], xa[3];
        for (int ix=0; ix<3; ++ix) { xi[ix]=coord[ix][i]; xa[ix]=cosurf[ix+1][ips]; }
        arat[i]+=cosurf[4][ips];
        double tm_flat[9];
        // Fortran tm(row=axis, col=component)锛汿mat 瀛樼殑鏄?[axis][component]锛?
        // mfinel 閲?TM(ix,j)=tm(ix,j)锛屾晠鎵撳寘鏃跺繀椤昏浆缃紙鍘熸潵鍐欐垚 [j2][ix] 鏄敊鐨勶級銆?
        for (int j2=0; j2<3; ++j2)
            for (int ix=0; ix<3; ++ix) tm_flat[j2*3+ix]=Tmat[i][ix][j2];
        int nfl1 = 0;
        mfinel(ips, 1, finel_flat.data(), nar_csm, nsetf, nset, rsc.data(),
               nipsrs, dirvec_flat.data(), tm_flat, xi, ri, nfl1, ioldcv,
               lenabc, maxrs);
        double aa = 0.0;
        for (int k=1; k<=nfl1; ++k) {
            double x1=finel_flat[(k-1)*4+0];
            double x2=finel_flat[(k-1)*4+1];
            double x3=finel_flat[(k-1)*4+2];
            double x4=finel_flat[(k-1)*4+3];
            aa += fdiagr * std::sqrt(x4*x4*x4);
            for (int l=1; l<=k-1; ++l) {
                double y1=finel_flat[(l-1)*4+0];
                double y2=finel_flat[(l-1)*4+1];
                double y3=finel_flat[(l-1)*4+2];
                double y4=finel_flat[(l-1)*4+3];
                aa += 2*x4*y4/std::sqrt((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2)+(x3-y3)*(x3-y3));
            }
        }
        amat[((ips+1)*ips)/2] = aa/(cosurf[4][ips]*cosurf[4][ips]);
        rdat[i] += aa;
        for (int jps=1; jps<=ips-1; ++jps) {
            int j=iatsp[jps];
            double d2 = 0.0;
            double xj[3] = {0.0, 0.0, 0.0};
            for (int ix=0; ix<3; ++ix) {
                xj[ix] = coord[ix][j];
                double xb = cosurf[ix+1][jps];
                d2 += (xb-xa[ix])*(xb-xa[ix]);
            }
            if (d2 > disex2) {
                aa = 1.0/std::sqrt(d2);
            } else {
                double rj = srad[j];
                double tm2_flat[9];
                for (int j2=0; j2<3; ++j2)
                    for (int ix=0; ix<3; ++ix) tm2_flat[j2*3+ix]=Tmat[j][ix][j2];
                int nfl2 = 0;
                mfinel(jps, 2, finel_flat.data(), nar_csm, nsetf, nset, rsc.data(),
                       nipsrs, dirvec_flat.data(), tm2_flat, xj, rj, nfl2, ioldcv,
                       lenabc, maxrs);
                aa = 0.0;
                for (int k=1; k<=nfl1; ++k) {
                    double x1=finel_flat[(k-1)*4+0];
                    double x2=finel_flat[(k-1)*4+1];
                    double x3=finel_flat[(k-1)*4+2];
                    double x4=finel_flat[(k-1)*4+3];
                    for (int l=1; l<=nfl2; ++l) {
                        double y1=finel_flat[1200+(l-1)*4+0];
                        double y2=finel_flat[1200+(l-1)*4+1];
                        double y3=finel_flat[1200+(l-1)*4+2];
                        double y4=finel_flat[1200+(l-1)*4+3];
                        aa += x4*y4/std::sqrt((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2)+(x3-y3)*(x3-y3));
                    }
                }
                aa = aa/(cosurf[4][ips]*cosurf[4][jps]);
            }
            amat[((ips-1)*ips)/2+jps] = aa;
            if (i == j) rdat[i] += 2*aa*cosurf[4][ips]*cosurf[4][jps];
        }
    }
    // PERFORM CHOLESKY FACTORIZATION (F90: call coscl1(amat, nsetf, nps, info);
    // nsetf is overwritten with the packed-offset table and reused by coscl2).
    if (!coscav_dump.empty()) {
        std::ofstream fp(coscav_dump + ".pre");
        int ntri_pre = nps * (nps + 1) / 2;
        for (int i = 1; i <= ntri_pre; ++i) {
            std::uint64_t b;
            std::memcpy(&b, &amat[i], 8);
            fp << std::hex << std::uppercase << std::setfill('0')
               << std::setw(16) << b << std::dec << '\n';
        }
    }
    {
        int info = 0;
        coscl1(amat.data(), nsetf, nps, info);
    }
    if (!coscav_dump.empty()) {
        std::ofstream f(coscav_dump + ".cpp");
        auto H = [&](double v) {
            std::uint64_t b;
            std::memcpy(&b, &v, 8);
            f << std::hex << std::uppercase << std::setfill('0')
              << std::setw(16) << b << std::dec << '\n';
        };
        f << nps << '\n';
        H(area); H(cosvol);
        int ntri = nps * (nps + 1) / 2;
        for (int i = 1; i <= ntri; ++i) H(amat[i]);
        for (int ix = 1; ix <= 4; ++ix)
            for (int ips = 1; ips <= nps; ++ips) H(cosurf[ix][ips]);
        for (int ips = 1; ips <= nps; ++ips) f << iatsp[ips] << '\n';
        for (int ips = 1; ips <= nps; ++ips) f << nar_csm[ips] << '\n';
        for (int ips = 1; ips <= nps; ++ips) f << nsetf[ips] << '\n';
        for (int a = 1; a <= numat; ++a) H(arat[a]);
    }
done_coscav:;
}

void mkbmat() {
    int iden = 0;
    for (int i = 1; i <= numat; ++i) {
        int ia = nfirst[i];
        int idel = nlast[i] - ia + 1;
        int nati = nat[i];
        double ddi = dd[nati] * a0;
        double qqi2 = (a0 * qq[nati]) * (a0 * qq[nati]);
        for (int ips = 1; ips <= nps; ++ips) {
            double xa1 = cosurf[1][ips]-coord[0][i];
            double xa2 = cosurf[2][ips]-coord[1][i];
            double xa3 = cosurf[3][ips]-coord[2][i];
            double dist = xa1*xa1+xa2*xa2+xa3*xa3;
            double rm1 = 1.0 / std::sqrt(dist);
            bmat[iden+1][ips] = rm1;
            if (idel > 1) {
                double rm3 = rm1*rm1*rm1;
                double rm5 = rm3*rm1*rm1;
                bmat[iden+3][ips]  = rm1 + 3*xa1*xa1*qqi2*rm5 - qqi2*rm3;
                bmat[iden+6][ips]  = rm1 + 3*xa2*xa2*qqi2*rm5 - qqi2*rm3;
                bmat[iden+10][ips] = rm1 + 3*xa3*xa3*qqi2*rm5 - qqi2*rm3;
                bmat[iden+2][ips] = xa1*ddi*rm3;
                bmat[iden+4][ips] = xa2*ddi*rm3;
                bmat[iden+7][ips] = xa3*ddi*rm3;
                bmat[iden+5][ips] = 3*xa1*xa2*qqi2*rm5;
                bmat[iden+8][ips] = 3*xa1*xa3*qqi2*rm5;
                bmat[iden+9][ips] = 3*xa3*xa2*qqi2*rm5;
                if (idel > 4) {
                    for (int iii = iden+11; iii <= iden+44; ++iii) bmat[iii][ips] = 0.0;
                    for (int iii = 5; iii <= 9; ++iii)
                        bmat[iden + (iii*(iii+1))/2][ips] = rm1;
                }
            }
        }
        iden += idel*(idel+1)/2;
    }

}

// Parse VDW(:El=r:...) keyword into per-element vdW radii; verify used radii.
void extvdw(double* vdw, const double* refvdw) {
    std::string line;
    size_t i = keywrd.find(" VDW(");
    if (i == std::string::npos) {
        line = " ";
    } else {
        if (keywrd[i + 5] != ';') keywrd.insert(i + 5, ";");
        size_t j = keywrd.find(")", i + 5);
        std::string seg = keywrd.substr(i + 5, j - (i + 5));
        for (char& ch : seg) if (ch == ':' || ch == ',') ch = ';';
        line = seg;
    }
    for (int t = 1; t <= 107; ++t) vdw[t] = refvdw[t];
    if (line != " ") {
        for (int el = 1; el <= 107; ++el) {
            int j = (elemts_C::cap_elemnt[el].size() >= 2 && elemts_C::cap_elemnt[el][1] == ' ') ? 1 : 2;
            std::string key = ";" + elemts_C::cap_elemnt[el].substr(0, j) + "=";
            size_t pos = line.find(key);
            if (pos != std::string::npos) vdw[el] = reada(line, (int)pos + 1);
        }
    }
    for (int at = 1; at <= numat; ++at) {
        if (nat[at] > 102) continue;
        if (vdw[nat[at]] > 900.0) {
            std::string msg = "MISSING VAN DER WAALS RADIUS " + elemts_C::elemnt[nat[at]];
            mopend(msg.c_str());
            return;
        }
    }
}

// Generate the near-regular spherical direction vectors (vertex tessellation).
// dirvec is (4, nppa) column-major: rows 1..3 = direction, row 4 = area weight.
void dvfill(int nppa, double* dirvec) {
    auto D=[&](int ix,int i) -> double& { return dirvec[(i-1)*4 + (ix-1)]; };
    auto Ds=[&](int ix,int i){ return &dirvec[(i-1)*4 + (ix-1)]; };
    int kset[31][3], fset[21][4];
    {
        int ks[60]={1,2,1,3,1,4,1,5,1,6,12,11,12,10,12,9,12,8,12,7,
                    2,3,3,4,4,5,5,6,6,2,7,8,8,9,9,10,10,11,11,7,2,
                    7,7,3,3,8,8,4,4,9,9,5,5,10,10,6,6,11,11,2};
        for(int i=0;i<60;++i){ int kk=i/2+1,kj=i%2; kset[kk][kj+1]=ks[i]; }
        int fs[60]={1,2,3,1,3,4,1,4,5,1,5,6,1,6,2,12,11,10,12,10,9,12,9,8,
                    12,8,7,12,7,11,2,3,7,3,4,8,4,5,9,5,6,10,6,2,11,7,8,3,
                    8,9,4,9,10,5,10,11,6,11,7,2};
        for(int i=0;i<60;++i){ int kk=i/3+1,kj=i%3; fset[kk][kj+1]=fs[i]; }
    }
    D(1,1)=-1; D(2,1)=0; D(3,1)=0;
    int nd=1;
    double r=std::sqrt(0.8), h=std::sqrt(0.2);
    for (int i=-1;i<=1;i+=2)
      for (int j=1;j<=5;++j) {
        ++nd;
        double beta = 1.0 + j*0.4*pi + (i+1)*0.1*pi;
        D(2,nd)=r*std::cos(beta); D(3,nd)=r*std::sin(beta); D(1,nd)=i*h;
      }
    D(1,12)=1; D(2,12)=0; D(3,12)=0;
    nd = 12;
    double cphi=std::cos(1.0), sphi=std::sin(1.0);
    for (int i=1;i<=12;++i) {
        double xx=D(1,i), yy=D(2,i);
        D(1,i)=cphi*xx+sphi*yy; D(2,i)=-sphi*xx+cphi*yy;
    }
    int m2=(nppa-2)/10;
    int m=(int)std::floor(std::sqrt((double)m2)+0.5);
    int k=0;
    if (m2 != m*m) { k=1; m2/=3; m=(int)std::floor(std::sqrt((double)m2)+0.5); }
    if (10*(k?3:1)*m*m+2 != nppa) return;
    for (int i=1;i<=30;++i) {
        int na=kset[i][1], nb=kset[i][2];
        for (int j=1;j<=m-1;++j) {
            ++nd;
            for (int ix=1;ix<=3;++ix)
                D(ix,nd)=D(ix,na)*(m-j)+D(ix,nb)*j;
        }
    }
    for (int i=1;i<=20;++i) {
        int na=fset[i][1], nb=fset[i][2], nc=fset[i][3];
        for (int j1=1;j1<=m-1;++j1)
          for (int j2=1;j2<=m-j1-1;++j2) {
            ++nd;
            for (int ix=1;ix<=3;++ix)
                D(ix,nd)=D(ix,na)*(m-j1-j2)+D(ix,nb)*j1+D(ix,nc)*j2;
          }
        if (k != 0) {
            double t=1.0/3;
            for (int j1=0;j1<=m-1;++j1)
              for (int j2=0;j2<=m-j1-1;++j2) {
                ++nd;
                for (int ix=1;ix<=3;++ix)
                    D(ix,nd)=D(ix,na)*(m-j1-j2-2*t)+D(ix,nb)*(j1+t)+D(ix,nc)*(j2+t);
              }
            t=2.0/3;
            for (int j1=0;j1<=m-2;++j1)
              for (int j2=0;j2<=m-j1-2;++j2) {
                ++nd;
                for (int ix=1;ix<=3;++ix)
                    D(ix,nd)=D(ix,na)*(m-j1-j2-2*t)+D(ix,nb)*(j1+t)+D(ix,nc)*(j2+t);
              }
        }
    }
    double sumar=0;
    for (int i=1;i<=nppa;++i) {
        double dist=0;
        for (int ix=1;ix<=3;++ix) dist += D(ix,i)*D(ix,i);
        dist = 1.0/std::sqrt(dist);
        double dist2 = (m*dist)*(m*dist);
        for (int ix=1;ix<=3;++ix) D(ix,i)*=dist;
        double ar = (i<=12) ? 5.0 : 6.0*dist2;
        D(4,i)=ar; sumar+=ar;
    }
    sumar = 4*pi/sumar;
    for (int i=1;i<=nppa;++i) D(4,i)*=sumar;
}

void addnuc() {
    double fcon = a0 * ev;
    for (int i = 1; i <= lm61; ++i) qdenet[i][1] = 0.0;
    for (int i = 1; i <= numat; ++i) qdenet[idenat[i]][1] = tore[nat[i]];
    for (int i = 1; i <= nps; ++i) {
        double phi = 0.0;
        for (int j = 1; j <= lm61; ++j) phi += bmat[j][i] * qdenet[j][1];
        phinet[i][1] = phi;
    }
    // coscl2 solves over 1..nps contiguous entries; the C++ storage of
    // qscnet/phinet is [segment][slot], so copy the slot-1 column to flat
    // 1-based buffers (F90 passes qscnet(1,1)/phinet(1,1) column-major).
    std::vector<double> q_flat(nps + 1, 0.0), p_flat(nps + 1, 0.0);
    for (int k = 1; k <= nps; ++k) p_flat[k] = phinet[k][1];
    std::vector<int> id_c2(nps + 1);
    id_c2[1] = 0;
    for (int i = 2; i <= nps; ++i) id_c2[i] = id_c2[i - 1] + i - 1;
    coscl2(amat.data(), id_c2, &q_flat[0], &p_flat[0], nps);
    double enclr = 0.0;
    for (int i = 1; i <= nps; ++i) {
        qscnet[i][1] = -fepsi * q_flat[i];
        enclr += qscnet[i][1] * phinet[i][1];
    }
    enuclr += fcon * enclr / 2.0;
}
int ijbo(int i, int j);

// COSMO initialization: radii, dielectric constants, segment counts, allocations.
void cosini(bool l_print) {
    static const double rvdw[108] = {
        0,
        1.30,1.64,2.13,2.19,2.05,2.00,1.83,1.72,1.72,1.80,2.66,2.02,2.41,2.46,2.11,2.16,
        2.05,2.20,3.22,2.54,2.64,2.64,2.52,2.40,2.46,2.41,2.40,1.91,1.64,1.63,2.19,2.46,
        2.22,2.22,2.16,2.36,3.78,3.44,3.39,3.33,3.28,2.57,2.57,2.57,2.57,1.91,2.01,1.85,
        2.26,2.54,2.53,2.41,2.32,2.53,4.00,3.47,2.81,2.81,2.81,2.81,2.81,2.81,2.81,2.81,
        2.81,2.81,2.81,2.81,2.81,2.81,2.81,2.57,2.57,2.57,2.57,2.57,2.57,2.05,1.94,1.81,
        2.29,2.36,2.64,2.64,2.63,2.69,2.57,2.57,2.57,2.57,2.57,2.18,2.57,2.57,2.57,2.57,
        2.57,2.57,2.57,2.57,2.57,2.57,2.0,2.0,2.0,2.0,2.0
    };
    std::vector<double> usevdw(108, 0.0);

    lenabc = std::max(100, nspa * numat);

    ipiden.assign(lm61 + 1, 0);
    idenat.assign(numat + 1, 0);
    gden.assign(lm61 + 1, 0.0);
    qdenet.assign(lm61 + 1, std::vector<double>(4, 0.0));
    phinet.assign(lenabc + 2, std::vector<double>(4, 0.0));
    qscnet.assign(lenabc + 2, std::vector<double>(4, 0.0));
    qscat.assign(numat + 1, 0.0);
    srad.assign(numat + 1, 0.0);
    nn.assign(4, std::vector<int>(numat + 1, 0));
    qden.assign(lm61 + 1, 0.0);
    iatsp.assign(lenabc + 2, 0);
    isude.assign(3, std::vector<int>(30 * numat + 1, 0));
    nar_csm.assign(lenabc + 2, 0);
    arat.assign(numat + 1, 0.0);
    sude.assign(3, std::vector<double>(30 * numat + 1, 0.0));
    nsetf.assign(lenabc + 2, 0);
    cosurf.assign(5, std::vector<double>(lenabc + 1, 0.0));
    if (!mozyme) {
        if (lenabc > 22000) {
            if (l_print) mopend("Data set too large to run using COSMO (try MOZYME).");
            moperr = true;
            return;
        }
        abcmat.assign(lenabc + 1, 0.0);
        xsp.assign(4, std::vector<double>(lenabc + 1, 0.0));
        nset.assign(std::max(1, nppa * numat) + 1, 0);
        bh.assign(lenabc + 1, 0.0);
        bmat.assign(lm61 + 1, std::vector<double>(lenabc + 1, 0.0));
        amat.assign((lenabc * (lenabc + 1)) / 2 + 1, 0.0);
        cmat.assign((lm61 * (lm61 + 1)) / 2 + 1, 0.0);
    }

    extvdw(usevdw.data(), rvdw);
    if (moperr) return;

    rsolv = 1.3;
    ioldcv = 0;
    int inrsol = (int)keywrd.find(" RSOLV=");
    if (inrsol >= 0) rsolv = reada(keywrd, inrsol);
    if (rsolv < 0.5) { /* RSOLV IS SET TO 0.5 */ return; }
    if (moperr) return;

    double ri1 = 2.0;
    int incif = (int)keywrd.find("N**2");
    if (incif >= 0) {
        ri1 = reada(keywrd, incif + 4);
        if (ri1 < 0.0) { mopend("N**2 CANNOT BE NEGATIVE"); return; }
    }
    fnsq = (ri1 - 1.0) / (ri1 + 0.5);
    nps = 0;
    rsolv = std::max(rsolv, 0.5);

    double disex = 4.0;
    int indise = (int)keywrd.find(" DISEX=");
    if (indise >= 0) disex = reada(keywrd, indise);

    int iden = 0;
    for (int i = 1; i <= numat; ++i) {
        int nfi = nfirst[i];
        idenat[i] = iden + 1;
        int idel = nlast[i] + 1 - nfi;
        if (mozyme) {
            int i0 = ijbo(i, i);
            for (int j = 1; j <= idel; ++j) {
                for (int k = 1; k <= j; ++k) {
                    ++iden; ++i0;
                    ipiden[iden] = i0;
                    gden[iden] = -2.0;
                }
                gden[iden] = -1.0;
            }
        } else {
            for (int j = 1; j <= idel; ++j) {
                int nfj = nfi - 1 + j;
                int i0 = (nfj * (nfj - 1)) / 2 + nfi - 1;
                for (int k = 1; k <= j; ++k) {
                    ++iden;
                    ipiden[iden] = i0 + k;
                    gden[iden] = -2.0;
                }
                gden[iden] = -1.0;
            }
        }
        srad[i] = usevdw[nat[i]];
    }

    n0[1] = nspa;
    if (nspa != 42) {
        double x = (nspa - 2) / 10.0 + 1e-8;
        int s1 = (int)std::sqrt(x), s2 = (int)std::sqrt(x / 3.0);
        int n1 = 10 * s1 * s1 + 2;
        int n2 = 30 * s2 * s2 + 2;
        n0[1] = std::max(std::max(n1, n2), 12);
    }
    nspa = n0[1];
    int n1b = (n0[1] - 2) / 10;
    n0[2] = 10 * n1b / 3 + 2;
    if (n1b % 3 != 0) { int s = (int)std::sqrt(n1b / 3.0); n0[2] = 10 * s * s + 2; }
    n0[2] = std::max(n0[2], 12);
    if (n0[1] + n0[2] > 1082) mopend("CHOSE NSPA < 0.75*1082");

    dvfill(n0[1], &dirsm[0][0]);
    // F90: dirsm(1, n0(1)+1) 鈥?column n0(1)+1 starts at planar offset n0(1)*4.
    dvfill(n0[2], &dirsm[0][0] + n0[1] * 4);
    disex2 = 4.0 * std::pow(1.7 * disex, 2.0) / nspa;
    dvfill(1082, &dirvec[0][0]);
}

// Close concave regions of the cavity by adding ring/intersection segments.
// surclo(): cosmo.F90:1202-1663 鐨勫畬鏁寸Щ妞嶏紙闂悎鍑归潰鍖哄煙鐨勭幆娈电敓鎴愶級銆?
// 绾﹀畾锛歝oord 骞抽摵涓?[ix][atom]锛坕x=0..2锛宎tom=1..numat锛夛紱rsc 涓哄垪涓诲簭
// 鎵佸钩缂撳啿锛屽厓绱?(ix,seg) 浣嶄簬 rsc[(seg-1)*4+(ix-1)]锛坕x=1..4锛夈€?
// cosurf/iatsp/nar_csm/nsetf/isude/sude/nps/lenabc/nipc 浣跨敤 cosmo_C 鍏ㄥ眬閲忋€?
void surclo(double* coord, const int* nipa, const int* lipa, const bool* din,
            int dim_din, double* rsc, int* isort, int* ipsrs, int* nipsrs,
            const int* nat, const double* srad, int maxrs) {
    (void)dim_din;
    const double twopi = 2.0 * pi;
    auto X = [&](int ix, int a) -> double { return coord[(ix - 1) * (numat + 1) + a]; };
    auto R = [&](int ix, int seg) -> double& { return rsc[(seg - 1) * 4 + (ix - 1)]; };
    auto C = [&](int ix, int ips) -> double& { return cosurf[ix][ips]; };

    // ---- differential-test dump (env SURCLO_DUMP=<prefix>) ----
    std::string surclo_dump;
    double surclo_area0 = area, surclo_cosvol0 = cosvol;
    if (const char* dp = std::getenv("SURCLO_DUMP")) {
        surclo_dump = dp;
        std::ofstream f(surclo_dump + ".in");
        f << std::setprecision(17);
        auto W = [&](double v) { f << v << '\n'; };
        W(numat); W(maxrs); W(lenabc); W(nps); W((double)n0[1]); W((double)n0[2]); W(rsolv);
        for (int a = 1; a <= numat; ++a)
            for (int ix = 1; ix <= 3; ++ix) W(coord[(ix - 1) * (numat + 1) + a]);
        for (int a = 1; a <= numat; ++a) W(srad[a]);
        for (int a = 1; a <= numat; ++a) W((double)nat[a]);
        for (int a = 1; a <= numat; ++a) W((double)nipa[a]);
        for (int a = 1; a <= numat; ++a) W(din[a] ? 1.0 : 0.0);
        for (int i = 1; i <= maxrs; ++i) W((double)lipa[i]);
        for (int ix = 1; ix <= 4; ++ix)
            for (int ips = 1; ips <= lenabc + 1; ++ips) W(C(ix, ips));
        for (int ips = 1; ips <= lenabc + 1; ++ips) W((double)iatsp[ips]);
        for (int ips = 1; ips <= lenabc + 1; ++ips) W((double)nar_csm[ips]);
        for (int ips = 1; ips <= lenabc + 1; ++ips) W((double)nsetf[ips]);
    }

    nipc = 0;
    int nrs = 0;
    for (int i = 1; i <= numat; ++i) nipsrs[i] = 0;
    int ilipa = 0;

    for (int ia = 1; ia <= numat - 1; ++ia) {
        if (!din[ia]) {
            double ra = srad[ia] + rsolv;
            double xta[4][4] = {};
            for (int ix = 1; ix <= 3; ++ix) xta[ix][1] = X(ix, ia);
            for (int iib = ilipa + 1; iib <= ilipa + nipa[ia]; ++iib) {
                int ib = lipa[iib];
                if (!(ib > ia)) continue;
                if (din[ib]) continue;
                double rb = srad[ib] + rsolv;
                double dab = 0.0;
                int nsab = 0;
                double xx[4] = {};
                for (int ix = 1; ix <= 3; ++ix) {
                    xta[ix][2] = X(ix, ib);
                    xx[ix] = xta[ix][2] - xta[ix][1];
                    dab += xx[ix] * xx[ix];
                }
                dab = std::sqrt(dab);
                double aa, ab, aar, abr, aad, abd, rinc;
                ansude(ra - rsolv, rb - rsolv, dab, rsolv,
                       aa, ab, aar, abr, aad, abd, rinc);
                double cosa = (ra * ra + dab * dab - rb * rb) / (2.0 * dab * ra);
                double cosb = (rb * rb + dab * dab - ra * ra) / (2.0 * dab * rb);
                double sina = std::sqrt(1.0 - cosa * cosa);
                double da = ra * cosa;
                double hh = ra * sina;
                double ddd = rsolv * (cosa + cosb) / dab;
                double fz[4] = {};
                fz[1] = (1.0 - std::cos(hh * pi / ra)) / 2.0;
                fz[2] = (1.0 - std::cos(hh * pi / rb)) / 2.0;
                if (cosa * cosb < 0.0) { fz[1] = 1.0; fz[2] = 1.0; }
                double yx[4] = {};
                yx[1] = rsolv / ra;
                yx[2] = rsolv / rb;
                double xd[4] = {};
                for (int ix = 1; ix <= 3; ++ix) xd[ix] = xta[ix][1] + da * xx[ix] / dab;
                double rvx[4] = {}, rvy[4] = {};
                rvx[1] = xx[2] * 3.0 - xx[3] * 2.0;
                rvx[2] = xx[3] * 1.0 - xx[1] * 3.0;
                rvx[3] = xx[1] * 2.0 - xx[2] * 1.0;
                double dist = std::sqrt(rvx[1] * rvx[1] + rvx[2] * rvx[2] + rvx[3] * rvx[3]);
                for (int ix = 1; ix <= 3; ++ix) rvx[ix] = hh * rvx[ix] / dist;
                rvy[1] = (xx[2] * rvx[3] - xx[3] * rvx[2]) / dab;
                rvy[2] = (xx[3] * rvx[1] - xx[1] * rvx[3]) / dab;
                rvy[3] = (xx[1] * rvx[2] - xx[2] * rvx[1]) / dab;

                int iset[51] = {};
                double phiset[51] = {};
                static double tarset[51] = {};  // F90 SAVE
                int ntrp = 0, ntrp2 = 0;

                for (int iic = ilipa + 1; iic <= ilipa + nipa[ia]; ++iic) {
                    int ic = lipa[iic];
                    if (ic == ib) continue;
                    double rc = srad[ic] + rsolv;
                    double dabc = 0.0, sp = 0.0;
                    for (int ix = 1; ix <= 3; ++ix) {
                        double xxx = X(ix, ic) - xd[ix];
                        xta[ix][3] = X(ix, ic);
                        sp += xxx * xx[ix];
                        dabc += xxx * xxx;
                    }
                    dabc = std::sqrt(dabc);
                    cosa = sp / dab / dabc;
                    sina = std::sqrt(std::max(1.0e-28, 1.0 - cosa * cosa));
                    double cj = (dabc * dabc + hh * hh - rc * rc) / (2.0 * dabc * hh * sina);
                    if (cj < 1.0) ntrp2 += 2;
                    if (cj <= 1.0 && cj >= -1.0) {
                        double sj = std::sqrt(1.0 - cj * cj);
                        double tvx[4] = {}, tvy[4] = {};
                        for (int ix = 1; ix <= 3; ++ix)
                            tvx[ix] = (xta[ix][3] - xd[ix]) - cosa * dabc * xx[ix] / dab;
                        dist = std::sqrt(tvx[1] * tvx[1] + tvx[2] * tvx[2] + tvx[3] * tvx[3]);
                        for (int ix = 1; ix <= 3; ++ix) tvx[ix] = hh * tvx[ix] / dist;
                        tvy[1] = (xx[2] * tvx[3] - xx[3] * tvx[2]) / dab;
                        tvy[2] = (xx[3] * tvx[1] - xx[1] * tvx[3]) / dab;
                        tvy[3] = (xx[1] * tvx[2] - xx[2] * tvx[1]) / dab;
                        for (int l = -1; l <= 1; l += 2) {
                            int il = ntrp + 1;
                            double trp[4] = {};
                            for (int ix = 1; ix <= 3; ++ix)
                                trp[ix] = xd[ix] + cj * tvx[ix] + sj * tvy[ix] * l;
                            bool skip = false;
                            for (int ik = ilipa + 1; ik <= ilipa + nipa[ia]; ++ik) {
                                int k = lipa[ik];
                                if (k != ib && k != ic) {
                                    double dabck = 0.0;
                                    for (int ix = 1; ix <= 3; ++ix) {
                                        double d = trp[ix] - X(ix, k);
                                        dabck += d * d;
                                    }
                                    dabck = std::sqrt(dabck);
                                    if (dabck < srad[k] + rsolv) { skip = true; break; }
                                }
                            }
                            if (skip) continue;
                            ++ntrp;
                            double spx = 0.0, spy = 0.0;
                            for (int ix = 1; ix <= 3; ++ix) {
                                spx += rvx[ix] * (trp[ix] - xd[ix]);
                                spy += rvy[ix] * (trp[ix] - xd[ix]);
                            }
                            double phi = std::acos(spx / (hh * hh + 1.0e-10));
                            if (spy < 0.0) phi = -phi;
                            phiset[il] = phi;
                            sp = 0.0;
                            for (int ix = 1; ix <= 3; ++ix)
                                sp += (-spy * rvx[ix] + spx * rvy[ix]) * (trp[ix] - xta[ix][3]);
                            iset[ntrp] = (sp < 0.0) ? -1 : 1;
                            double ee[4][4] = {};
                            sp = 0.0;
                            double dac = 0.0, dbc = 0.0;
                            for (int ix = 1; ix <= 3; ++ix) {
                                ee[ix][1] = trp[ix] + rsolv / srad[ia] * (xta[ix][1] - trp[ix]);
                                ee[ix][2] = trp[ix] + rsolv / srad[ib] * (xta[ix][2] - trp[ix]);
                                ee[ix][3] = trp[ix] + rsolv / srad[ic] * (xta[ix][3] - trp[ix]);
                                sp += (ee[ix][1] - ee[ix][3]) * (ee[ix][2] - ee[ix][3]);
                                dac += (ee[ix][1] - ee[ix][3]) * (ee[ix][1] - ee[ix][3]);
                                dbc += (ee[ix][2] - ee[ix][3]) * (ee[ix][2] - ee[ix][3]);
                            }
                            tarset[il] = 0.8 * std::sqrt(dac * dbc - sp * sp) / 12.0;
                        }
                    }
                }

                if ((ntrp % 2) != 0) mopend("ODD NTRP");
                if (ntrp > 18) mopend("NTRP TOO LARGE");
                if (ntrp + ntrp2 == 0) {
                    phiset[1] = 0.0;
                    phiset[2] = twopi;
                    tarset[2] = 0.0;
                    iset[1] = 1;
                    iset[2] = -1;
                    ntrp = 2;
                }
                for (;;) {
                    int icnt = 0;
                    for (int l = 2; l <= ntrp; ++l) {
                        if (phiset[l] < phiset[l - 1]) {
                            std::swap(phiset[l], phiset[l - 1]);
                            std::swap(iset[l], iset[l - 1]);
                            std::swap(tarset[l], tarset[l - 1]);
                            ++icnt;
                        }
                    }
                    if (icnt <= 0) {
                        if (iset[1] == -1) phiset[1] += twopi;
                        else break;
                    }
                }

                double sumphi = 0.0;
                for (int l = 2; l <= ntrp; l += 2) {
                    int k = l - 1;
                    double phiu = phiset[k];
                    double phio = phiset[l];
                    int nsa = (int)((phio - phiu) / 2.0 / pi * 20.0);
                    nsa = std::max(nsa + 1, 2);
                    sumphi += phio - phiu;
                    double dp = (phio - phiu) / (nsa - 1);
                    for (int ich = 1; ich <= 2; ++ich) {
                        int iat = (ich == 2) ? ia : ib;
                        double htr = (iat == ia) ? aar / twopi : abr / twopi;
                        for (int ja = ich; ja <= nsa; ja += 2) {
                            int jb = std::max(ja - 1, 1);
                            int jc = std::min(ja + 1, nsa);
                            double cphi, sphi, phi;
                            double xja[4] = {}, xjb[4] = {}, xjc[4] = {};
                            phi = phiu + (ja - 1) * dp;
                            cphi = std::cos(phi);
                            sphi = std::sin(phi);
                            for (int ix = 1; ix <= 3; ++ix) {
                                double ca = xd[ix] + (cphi * rvx[ix] + sphi * rvy[ix]) * fz[ich];
                                ca = ca + (xta[ix][ich] - ca) * yx[ich];
                                xja[ix] = ca + (xta[ix][3 - ich] - xta[ix][ich]) * ddd * (1.0 - fz[ich]);
                            }
                            phi = phiu + (jb - 1) * dp;
                            cphi = std::cos(phi);
                            sphi = std::sin(phi);
                            for (int ix = 1; ix <= 3; ++ix) {
                                double ca = xd[ix] + cphi * rvx[ix] + sphi * rvy[ix];
                                xjb[ix] = ca + (xta[ix][3 - ich] - ca) * yx[3 - ich];
                            }
                            phi = phiu + (jc - 1) * dp;
                            cphi = std::cos(phi);
                            sphi = std::sin(phi);
                            for (int ix = 1; ix <= 3; ++ix) {
                                double ca = xd[ix] + cphi * rvx[ix] + sphi * rvy[ix];
                                xjc[ix] = ca + (xta[ix][3 - ich] - ca) * yx[3 - ich];
                            }
                            ++nrs;
                            ++nsab;
                            double sp2 = 0.0, dist2 = 0.0, spn = 0.0, spn2 = 0.0;
                            ipsrs[nrs] = (ich == 2) ? ia : ib;
                            int iat2 = ipsrs[nrs];
                            ++nipsrs[iat2];
                            for (int ix = 1; ix <= 3; ++ix) {
                                R(ix, nrs) = (xja[ix] * 0.5 + xjb[ix] + xjc[ix]) / 2.5;
                                int i2 = (ix % 3) + 1;
                                int i3 = (i2 % 3) + 1;
                                double cnrs = (xjc[i2] - xjb[i2]) * (xja[i3] - xjb[i3]) -
                                              (xja[i2] - xjb[i2]) * (xjc[i3] - xjb[i3]);
                                dist2 += cnrs * cnrs;
                                spn += cnrs * (R(ix, nrs) - X(ix, iat2));
                                spn2 += cnrs * R(ix, nrs);
                            }
                            dist2 = 1.0 / std::sqrt(dist2);
                            if (spn < 0.0) dist2 = -dist2;
                            R(4, nrs) = (jc - jb) * dp * htr;
                            if (ja == 1) R(4, nrs) += tarset[k] * rinc;
                            if (ja == nsa) R(4, nrs) += tarset[l] * rinc;
                            cosvol += R(4, nrs) * spn2 * dist2;
                            area += R(4, nrs);
                        }
                    }
                }
                if (sumphi > 1.0e-10) {
                    ++nipc;
                    isude[1][nipc] = ia;
                    isude[2][nipc] = ib;
                    sumphi = sumphi / twopi;
                    sude[1][nipc] = aad * sumphi;
                    sude[2][nipc] = abd * sumphi;
                }
            }
        }
        ilipa += nipa[ia];
    }
    if (nrs > maxrs) mopend("NRS .GT. MAXRS IN SURCLO");

    // SORT THE SEGMENTS WITH RESPECT TO ATOMS
    int isum = 0;
    for (int iat = 1; iat <= numat; ++iat) {
        int isum2 = isum + nipsrs[iat];
        nipsrs[iat] = isum;
        isum = isum2;
    }
    for (int i = 1; i <= nrs; ++i) {
        int iat = ipsrs[i];
        int isum2 = nipsrs[iat] + 1;
        nipsrs[iat] = isum2;
        isort[i] = isum2;
    }
    for (int i = 1; i <= nrs; ++i) {
        while (isort[i] != i) {
            int is = isort[i];
            for (int ix = 1; ix <= 4; ++ix) std::swap(R(ix, i), R(ix, is));
            std::swap(ipsrs[i], ipsrs[is]);
            isort[i] = isort[is];
            isort[is] = is;
        }
    }

    // FIND FOR EACH RING-SEGMENT THE NEAREST SEGMENT OF THE SAME ATOM
    for (int i = 1; i <= nps; ++i) nipsrs[i] = 0;
    int iat0 = 0, ips1 = 0, ips0 = 1, ipsmin = 0;
    double d2max = 0.0;
    for (int i = 1; i <= nrs; ++i) {
        int iat = ipsrs[i];
        if (iat > iat0) {
            iat0 = iat;
            d2max = 16.0 * srad[iat] * srad[iat] / n0[1];
            if (nat[iat] == 1) d2max = d2max * n0[1] / n0[2];
            int ips = ips1 + 1;
            for (; ips <= nps; ++ips) if (iatsp[ips] == iat) break;
            ips0 = ips;
            for (ips = ips0; ips <= nps - 1; ++ips) if (iatsp[ips + 1] > iat) break;
            ips1 = ips;
        }
        double d2min = 1.0e6;
        for (int ips = ips0; ips <= ips1; ++ips) {
            double d2 = 0.0;
            for (int ix = 1; ix <= 3; ++ix) {
                double d = C(ix, ips) - R(ix, i);
                d2 += d * d;
            }
            if (d2 < d2min) { d2min = d2; ipsmin = ips; }
        }
        if (d2min > d2max) {
            ++ips1;
            for (int ip1 = nps; ip1 >= ips1; --ip1) {
                int ip2 = ip1 + 1;
                for (int ix = 1; ix <= 4; ++ix) C(ix, ip2) = C(ix, ip1);
                iatsp[ip2] = iatsp[ip1];
                nsetf[ip2] = nsetf[ip1];
                nar_csm[ip2] = nar_csm[ip1];
                nipsrs[ip2] = nipsrs[ip1];
            }
            ipsrs[i] = ips1;
            ++nps;
            for (int ix = 1; ix <= 4; ++ix) C(ix, ips1) = R(ix, i);
            iatsp[ips1] = iatsp[ips1 - 1];
            nar_csm[ips1] = 0;
            nipsrs[ips1] = 1;
            nsetf[ips1] = nsetf[ips1 + 1];
        } else {
            ipsrs[i] = ipsmin;
            ++nipsrs[ipsmin];
            double arseg = C(4, ipsmin);
            double arsegn = arseg + R(4, i);
            for (int ix = 1; ix <= 3; ++ix)
                C(ix, ipsmin) = (arseg * C(ix, ipsmin) + R(4, i) * R(ix, i)) / arsegn;
            C(4, ipsmin) = arsegn;
        }
    }

    // SORT THE RING-SEGMENTS WITH RESPECT TO PRIMARY SEGMENTS
    isum = 0;
    for (int ips = 1; ips <= nps; ++ips) {
        int isum2 = isum + nipsrs[ips];
        nipsrs[ips] = isum;
        isum = isum2;
    }
    for (int i = 1; i <= nrs; ++i) {
        int ips = ipsrs[i];
        int isum2 = nipsrs[ips] + 1;
        nipsrs[ips] = isum2;
        isort[i] = isum2;
    }
    for (int i = 1; i <= nrs; ++i) {
        while (isort[i] != i) {
            int is = isort[i];
            for (int ix = 1; ix <= 4; ++ix) std::swap(R(ix, i), R(ix, is));
            isort[i] = isort[is];
            isort[is] = is;
        }
    }

    if (!surclo_dump.empty()) {
        std::ofstream f(surclo_dump + ".cpp");
        auto H = [&](double v) {
            std::uint64_t b;
            std::memcpy(&b, &v, 8);
            f << std::hex << std::uppercase << std::setfill('0')
              << std::setw(16) << b << std::dec << '\n';
        };
        f << nps << '\n' << nipc << '\n' << nrs << '\n';
        H(area - surclo_area0);
        H(cosvol - surclo_cosvol0);
        for (int seg = 1; seg <= nrs; ++seg)
            for (int ix = 1; ix <= 4; ++ix) H(R(ix, seg));
        for (int i = 1; i <= nps; ++i) f << nipsrs[i] << '\n';
        for (int i = 1; i <= nrs; ++i) f << ipsrs[i] << '\n';
        for (int i = 1; i <= nrs; ++i) f << isort[i] << '\n';
        for (int ix = 1; ix <= 4; ++ix)
            for (int ips = 1; ips <= nps; ++ips) H(C(ix, ips));
        for (int ips = 1; ips <= nps; ++ips) f << iatsp[ips] << '\n';
        for (int ips = 1; ips <= nps; ++ips) f << nar_csm[ips] << '\n';
        for (int ips = 1; ips <= nps; ++ips) f << nsetf[ips] << '\n';
    }
}




