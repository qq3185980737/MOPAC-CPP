// ============================================================================
// 【中文导读】cosmo_C.cpp —— COSMO 模块全局量的存储定义
// ----------------------------------------------------------------------------
// 对应 cosmo_C.F90 的 module 数据：nps/lenabc/nipa 等计数、srad/cosurf/amat/
// bmat/phinet/qscnet 等数组、fepsi/ediel/area/cosvol 等标量。
//   - n0[] 为 1-based（n0[1]=nspa，n0[2]=氢方向数），代码里不要写 n0[0]；
//   - dirvec/dirsm 为 Fortran 列主序的 4 分量方向表（第 4 分量为权重）；
//   - amat 为打包下三角、cosurf(4,lenabc+1) 第 4 行存段面积。
// 这里的定义必须与 cosmo.cpp 的 1-based 索引约定严格一致。
// ============================================================================

// ============================================================================
// 【中文导读】cosmo_C.cpp —— COSMO 模块全局量的存储定义
// ----------------------------------------------------------------------------
// 对应 cosmo_C.F90 的 module 数据：nps/lenabc/nipa 等计数、srad/cosurf/amat/
// bmat/phinet/qscnet 等数组、fepsi/ediel/area/cosvol 等标量。
//   - n0[] 为 1-based（n0[1]=nspa，n0[2]=氢方向数），代码里不要写 n0[0]；
//   - dirvec/dirsm 为 Fortran 列主序的 4 分量方向表（第 4 分量为权重）；
//   - amat 为打包下三角、cosurf(4,lenabc+1) 第 4 行存段面积。
// 这里的定义必须与 cosmo.cpp 的 1-based 索引约定严格一致。
// ============================================================================

// cosmo_C.cpp — C++ translation of MOPAC 2016 "cosmo_C.F90".
#include "cosmo_C.h"

namespace cosmo_C {
bool iseps = false, noeps = false, useps = false, lpka = false;
int nspa = 0, nps = 0, nps2 = 0, nden = 0, lenabc = 0, nppa = 1082;
int amat_dim = 0, isude_dim = 0, nipc = 0, ioldcv = 0;
int n0[3] = {0, 0, 0};
std::vector<int> iatsp, nar_csm, nsetf, cosmo_i, ipiden, idenat, nset;
std::vector<std::vector<int>> isude, nn;
double fepsi = 0.0, rds = 0.0, disex2 = 0.0, ediel = 0.0, solv_energy = 0.0;
double area = 0.0, cosvol = 0.0;
double cif1 = 0.0, cif2 = 0.0, cdiagi = 0.0, fnsq = 0.0, rsolv = 0.0;
double tm[4][4][1083] = {};
double dirsm[5][1083] = {}, dirvec[5][1083] = {};
std::vector<double> amat, cmat, gden, qscat, arat, srad, abcmat, qden, bh, cdiag;
std::vector<std::vector<double>> bmat, phinet, qscnet, qdenet, cosurf, sude, xsp, cxy;
}
