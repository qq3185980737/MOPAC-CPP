// ============================================================================
// 【中文导读】esp_C.cpp —— 偶极矩/极化率/静电与介电
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 esp_C.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// esp_C.cpp — C++ translation of MOPAC 2016 "esp_C.F90" (ESP fit module data).
#include "esp_C.h"

namespace esp_C {
int nesp, idip, iz, ipx, isc, is_esp, icd, ipe, npr, ic, ip, ncc;
double dens, scale, cf, rms, rrms, dx, dy, dz, den;

double fv[9][822] = {};
double dex[98] = {};       // dex(-1:96)
double tf[3] = {33.0, 37.0, 41.0};
double fac[8] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0};

// S    P             D
// s  x y z   xx yy zz xy xz yz
int ixn[11] = {0, 0, 4, 0, 0, 8, 0, 0, 4, 4, 0};
int iyn[11] = {0, 0, 0, 4, 0, 0, 8, 0, 4, 0, 4};
int izn[11] = {0, 0, 0, 0, 4, 0, 0, 8, 0, 4, 4};
int jxn[11] = {0, 0, 1, 0, 0, 2, 0, 0, 1, 1, 0};
int jyn[11] = {0, 0, 0, 1, 0, 0, 2, 0, 1, 0, 1};
int jzn[11] = {0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 1};

std::vector<int> ind, itemp, ird, indc;
std::vector<double> cc, ex, temp, rad, es, b_esp, esp_array, cesp, cespml,
    al, qesp, td, dx_array, dy_array, dz_array, ptd, pexs, pce, pexpn, pewcx,
    pewcy, pewcz, pf0, pf1, pf2, pmtd, fc, qsc, espx, espp;
std::vector<std::vector<int>> iam;
std::vector<std::vector<double>> cespm, cen, potpc, co, potpt, ovl, cespm2,
    a2, exs, ce, expn, ewcx, ewcy, ewcz, f0, f1, u_esp, rnai, rnai1, rnai2,
    espi, exsr;
std::vector<std::vector<bool>> cequiv;
}  // namespace esp_C
