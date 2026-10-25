// ============================================================================
// 【中文导读】cosmo.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// cosmo.h — C++ translation of MOPAC 2016 "cosmo.F90" (small routines; heavy geometric ones stubbed).
#pragma once
#include <vector>

// Analytic area of two intersecting spheres (radii ra, rb, distance d, probe rs).
void ansude(double ra, double rb, double d, double rs,
            double& aar, double& abr, double& ara, double& arad,
            double& arb, double& arbd, double& rinc);

void ciint(const double* c34, double* pq34);
void addfck(double* pin, double* fin);
void addhcr();
void addnuc();

// Packed lower-triangle Cholesky factorization / solve.
void coscl1(double* a, std::vector<int>& id, int n, int& info);
void coscl2(const double* a, const std::vector<int>& id,
            double* x, const double* y, int n);
void diegrd(double* dxyz);
void dmecip(double* coeffs, double* deltap, double* delta, double* eig,
            double* vectci, const double* conf);
void mkbmat();
void extvdw(double* vdw, const double* refvdw);
void dvfill(int nppa, double* dirvec);
void mfinel(int ips, int k, double* finel,
            const std::vector<int>& nar_csm, const std::vector<int>& nsetf,
            const std::vector<int>& nset, double* rsc,
            const std::vector<int>& nipsrs,
            double* dirvec, double* tm, const double* x, double r,
            int& nfl, int ioldcv, int lenabc, int maxrs);
void coscav();
void cosini(bool l_print);
void surclo(double* coord, const int* nipa, const int* lipa, const bool* din,
            int dim_din, double* rsc, int* isort, int* ipsrs, int* nipsrs,
            const int* nat, const double* srad, int maxrs);
