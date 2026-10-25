// ============================================================================
// 【中文导读】mult33.cpp —— 数值/线性代数与几何工具
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mult33.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mult33.cpp
#include "mult33.h"
#include "symmetry_C.h"

using namespace symmetry_C;

void mult33(const double* fmat, int iplace) {
    // fmat is Fortran 3x3 column-major: fmat(i,l) = fmat[(l-1)*3 + (i-1)].
    double help[4][4] = {};
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j) {
            help[i][j] = 0.0;
            for (int k = 1; k <= 3; ++k) {
                double s = 0.0;
                for (int l = 1; l <= 3; ++l)
                    s += fmat[(l-1)*3 + (i-1)] * fmat[(k-1)*3 + (j-1)] * elem[l][k][iplace];
                help[i][j] += s;
            }
        }
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j)
            elem[i][j][iplace] = help[i][j];
}
