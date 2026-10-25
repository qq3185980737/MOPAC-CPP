// ============================================================================
// 【中文导读】write_unit_cell_hof.cpp —— 周期体系单胞生成热输出
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 moldat.F90:1118-1173 的 write_unit_cell_HOF 子程序。
// 作用：把总生成热换算成“每个单胞的生成热”，并报告参与统计的单胞数 z。
//   z 的决定方式（与 Fortran 完全一致）：
//     1) 关键字给了 Z=：z = mers(1)*mers(2)*mers(3)*Z，BCC 时再除 2；
//     2) 否则由当前分子的经验式求出最小公倍数式单胞数：
//        统计各元素原子数，求最小原子数 k，取最小 i 使 i*nel/k 全为整数，
//        单胞数 z = k/i。
// 输出格式：10X,'H.o.F. per unit cell    =',F17.5,' KCAL, for',I4,' unit cells'
// 仅在 id==3（周期体系）时由 writmo 调用。
// ============================================================================

// write_unit_cell_hof.cpp — C++ translation of moldat.F90:1118-1173
// ("write_unit_cell_HOF").
#include "write_unit_cell_hof.h"
#include "chanel_C.h"

#include <cmath>
#include <cstdio>
#include <string>

#include "common_arrays_C.h"
#include "molkst_C.h"
#include "reada.h"
#include "to_screen.h"

using namespace common_arrays_C;
using namespace molkst_C;

void write_unit_cell_HOF(int iprt) {
    if (iprt < 0) return;

    // ---- 1) 决定单胞数 z ----
    int z = 0;
    const std::size_t iz = keywrd.find(" Z=");
    if (iz != std::string::npos) {
        const int zv = static_cast<int>(std::lround(reada(keywrd, static_cast<int>(iz))));
        z = mers[1] * mers[2] * mers[3] * zv;
        if (keywrd.find(" BCC") != std::string::npos) z = z / 2;
    } else {
        int nel[108] = {0};
        for (int i = 1; i <= numat; ++i) nel[nat[i]] = nel[nat[i]] + 1;
        int j = 0;
        for (int i = 1; i <= 107; ++i) {
            if (nel[i] > 0) {
                ++j;
                nel[j] = nel[i];   // 压缩：前 j 个位置只留非零元素计数
            }
        }
        int k = 1000;
        for (int i = 1; i <= j; ++i)
            if (nel[i] < k) k = nel[i];
        // k 是经验式中最小原子数；找最小 i 使所有 i*nel/k 为整数
        int m = 0, i = 1;
        for (i = 1; i <= 10; ++i) {
            m = 0;
            for (int l = 1; l <= j; ++l)
                if (std::fabs((i * nel[l]) / static_cast<double>(k) -
                              (i * 1.0 * nel[l]) / k) > 1.0e-5)
                    m = 1;
            if (m == 0) break;
        }
        z = k / i;   // 经验单元数
    }
    if (z <= 0) return;   // 防御：避免除零

    // ---- 2) 组装并输出 ----
    char buf[256];
    std::snprintf(buf, sizeof(buf),
                  "          H.o.F. per unit cell    =%17.5f KCAL, for%4d unit cells",
                  escf / z, z);
    line = buf;
    if (iprt == 0) {
        to_screen(line);
    } else {
        std::fprintf(chanel_C::out_unit(), "%s\n", line.c_str());
    }
}
