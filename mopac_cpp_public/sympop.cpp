// ============================================================================
// 【中文导读】sympop.cpp —— 对称性/特征标
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 sympop.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// sympop.cpp — C++ translation of MOPAC 2016 "sympop.F90".
// Find a symmetry-related atom (via ipo) that is suitable for a transition
// dipole calculation and symmetrize the Hessian block for atom i.
#include "sympop.h"
#include "symmetry_C.h"
#include "symh.h"

using symmetry_C::nsym;
using symmetry_C::ipo;

void sympop(double* h, int i, int& iskip, double* deldip) {
    for (int j = 1; j <= nsym; ++j) {
        if (ipo[i][j] >= i) continue;
        symh(h, deldip, i, j, nullptr);
        iskip = 3;  // atom ipo(i,j) is suitable for transition dipole calc'n
        return;
    }
    iskip = 0;
}
