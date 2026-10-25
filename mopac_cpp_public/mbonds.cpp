// ============================================================================
// 【中文导读】mbonds.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mbonds.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mbonds.cpp — C++ translation of MOPAC 2016 "mbonds.F90".
#include "mbonds.h"
#include <algorithm>
#include <cmath>
#include "chanel_C.h"
#include "molkst_C.h"
int ijbo(int i, int j);
void mopend(const char* msg);
void mbonds(int locc, int lvir, const double* f, double* catom,
            const int* nfirst, int ii, int jj, bool* ui, bool* uj,
            bool& lok, const int* iorbs, double* cocc, double* cvir,
            int cocc_dim, int cvir_dim, int numat, int norbs, int morb,
            int mpack1) {
    (void)mpack1;
    lok = false;
    int iorbsi = iorbs[ii];
    int iorbsj = iorbs[jj];
    int ni = nfirst[ii] - 1;
    int nj = nfirst[jj] - 1;
    double summin = 0.0, summax = -1.0e5;
    bool ok = false;
    int i1 = 0, i2 = 0, j1 = 0, j2 = 0;
    for (int i = 1; i <= iorbsi; ++i) {
        if (!ui[i]) {
            for (int j = 1; j <= iorbsj; ++j) {
                if (!uj[j]) {
                    ok = true;
                    double sum = 0.0;
                    int l = ijbo(ii, jj);
                    for (int m = 1; m <= iorbsj; ++m)
                        for (int k = 1; k <= iorbsi; ++k) {
                            l = l + 1;
                            sum += catom[(i + ni - 1) * morb + (k - 1)] * f[l] *
                                   catom[(j + nj - 1) * morb + (m - 1)];
                        }
                    if (sum < summin) { i1=i; j1=j; summin=sum; }
                    if (sum > summax) { i2=i; j2=j; summax=sum; }
                }
            }
        }
    }
    if (!ok) return;
    lok = true;
    int k = 0;
    for (int i=1;i<=iorbsi;++i) if (ui[i]) ++k;
    for (int i=1;i<=iorbsj;++i) if (uj[i]) ++k;
    if ((locc+iorbsi+iorbsj-k)>cocc_dim || (lvir+iorbsi+iorbsj-k)>cvir_dim) {
        mopend("VALUE OF NLMO IS TOO SMALL");
        return;
    }
    double one, e12;
    if (summax > -summin) { i1=i2; j1=j2; one=-1.0; e12=-summax; }
    else { one=1.0; e12=summin; }
    double e11=0.0, e111=0.0;
    int l = ijbo(ii, ii);
    for (int i=1;i<=iorbsi;++i) {
        for (int j=1;j<=i-1;++j) {
            ++l;
            e111 += catom[(ni+i1-1)*morb+(i-1)]*f[l]*catom[(ni+i1-1)*morb+(j-1)];
        }
        ++l;
        e11 += catom[(ni+i1-1)*morb+(i-1)]*f[l]*catom[(ni+i1-1)*morb+(i-1)];
    }
    e11 += e111*2.0;
    l = ijbo(jj, jj);
    double e22=0.0, e221=0.0;
    for (int i=1;i<=iorbsj;++i) {
        for (int j=1;j<=i-1;++j) {
            ++l;
            e221 += catom[(nj+j1-1)*morb+(i-1)]*f[l]*catom[(nj+j1-1)*morb+(j-1)];
        }
        ++l;
        e22 += catom[(nj+j1-1)*morb+(i-1)]*f[l]*catom[(nj+j1-1)*morb+(i-1)];
    }
    e22 += e221*2.0;
    double d = e11-e22;
    double e = std::sqrt(4.0*e12*e12+d*d);
    if (d<0) e=-e;
    double alpha = std::min(0.866, std::max(std::sqrt(0.5*(1.0+d/e)), 0.5));
    double beta = -std::sqrt(1.0-alpha*alpha);
    if (e12<0) beta=-beta;
    for (int i=1;i<=iorbsi;++i) {
        cocc[locc+i] = catom[(ni+i1-1)*morb+(i-1)]*alpha;
        cvir[lvir+i] = catom[(ni+i1-1)*morb+(i-1)]*beta;
    }
    int lim = std::min(iorbsj, std::min(cocc_dim-locc-iorbsi, cvir_dim-lvir-iorbsi));
    for (int i=1;i<=lim;++i) {
        cocc[locc+iorbsi+i] = catom[(nj+j1-1)*morb+(i-1)]*beta*one;
        cvir[lvir+iorbsi+i] = -catom[(nj+j1-1)*morb+(i-1)]*alpha*one;
    }
    ui[i1] = true;
    uj[j1] = true;
}