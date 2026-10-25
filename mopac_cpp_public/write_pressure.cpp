// ============================================================================
// 【中文导读】write_pressure.cpp —— 周期体系“约束平移向量所需压力”
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 moldat.F90:1174-1248 的 write_pressure 子程序。
// 原理：对每个 Tv 原子（labels==107）收集它的 3 个平移向量分量上的梯度
// dsum 与分量值 dsum1，把梯度投影到 Tv 方向：
//     xi = (dsum · dsum1) / |dsum1|
// 再换算成压力（GPa）：
//     xi = -(4184*10^30)/fpc_10 * xi*|dsum1| / V
// 最后减去外压项并按 1e-9 缩放写入 press(ndim)，输出：
//     10X,"Tv(",I4,")  Pressure:",F7.2," GPa"
// 若使用内坐标（na(k)!=0），无法计算，打印两条提示并返回。
// ============================================================================

// write_pressure.cpp — C++ translation of moldat.F90:1174-1248
// ("write_pressure").
#include "write_pressure.h"
#include "chanel_C.h"

#include <cmath>
#include <cstdio>
#include <string>

#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "to_screen.h"
#include "volume.h"

using namespace common_arrays_C;
using namespace funcon_C;
using namespace molkst_C;

namespace {

// 3 维点积（对应 F90 ddot(3,...)）。
double dot3(const double* a, const double* b) {
    return a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

// 输出一行：iprt==0 走 to_screen（无前导空格）；iprt>0 走 list-directed
// 写文件（Fortran 会在字符值前补一个空格）。
void emit_line(int iprt, const std::string& s) {
    if (iprt == 0) {
        to_screen(s);
    } else {
        std::fprintf(chanel_C::out_unit(), " %s\n", s.c_str());
    }
}

}  // namespace

void write_pressure(int iprt) {
    if (iprt < 0) return;

    int m = 0, i1 = 0, i2 = 0, ndim = 0;
    double dsum[4] = {0.0, 0.0, 0.0, 0.0};
    double dsum1[4] = {0.0, 0.0, 0.0, 0.0};

    for (int i = 1; i <= nvar; ++i) {
        const int k = loc[1][i];
        const int l = labels[k];
        double xi = xparam[i];
        if (l == 107 && (m == 0 || k == i1)) {
            // 该变量属于 Tv 原子：累积 3 个平移分量
            i1 = k;
            ++m;
            dsum[m] = grad[i];
            dsum1[m] = xparam[i];
            if (m == 3) {
                if (na[k] != 0) {
                    // 内坐标下无法把梯度投影到平移向量
                    if (iprt == 0) {
                        to_screen("The pressure required to constrain translation vectors");
                        to_screen("can only be calculated if Cartesian coordinates are used.");
                    } else {
                        if (std::fabs(pressure) > 0.01) {
                            std::fprintf(chanel_C::out_unit(), " The pressure required to constrain translation vectors\n");
                            std::fprintf(chanel_C::out_unit(), " can only be calculated if Cartesian coordinates are used.\n");
                        }
                    }
                    return;
                }
                // 梯度在 Tv 方向上的标量分量
                xi = dot3(dsum, dsum1) / std::sqrt(dot3(dsum1, dsum1));
                // 力/面积 → 压力；体积 V 来自三个晶格矢量
                double vecs[9];
                for (int ii = 1; ii <= 3; ++ii)
                    for (int jj = 1; jj <= 3; ++jj)
                        vecs[(jj - 1) * 3 + (ii - 1)] = tvec[ii][jj];
                const double vol = volume(vecs, 3);
                xi = -(4184.0 * 1.0e30) / fpc_10 * xi *
                     std::sqrt(dot3(dsum1, dsum1)) / vol;
                if (std::fabs(xi) < 1.0e-20) continue;  // 梯度为零时不打印
                if (i2 == 0) {
                    emit_line(iprt,
                              "          Pressure required to constrain translation vectors");
                    i2 = 1;
                }
                xi = (xi - pressure * (4184.0 * 1.0e30) / fpc_10) * 1.0e-9;
                ++ndim;
                if (ndim <= 3) press[ndim] = xi;   // F90 press(3)
                char buf[128];
                std::snprintf(buf, sizeof(buf), "          Tv(%4d)  Pressure:%7.2f GPa", k, xi);
                emit_line(iprt, buf);
                m = 0;
            }
        }
    }
}
