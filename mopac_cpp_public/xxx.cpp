// ============================================================================
// 【中文导读】xxx.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 xxx.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// xxx.cpp
#include "xxx.h"
void xxx(char type, int i, int j, int k, int l, std::string& r) {
    // Fortran r*13: r(1)=type, digits follow. C++ 0-based: r[0]=type.
    r.clear();
    r.push_back(type);
    int ijk[4] = {i, j, k, l};
    for (int loop = 0; loop < 4; ++loop) {
        int ii = ijk[loop];
        if (ii == 0) continue;
        int i2 = ii / 10;
        if (i2 != 0) { r.push_back(char('0' + i2)); ii -= i2 * 10; }
        r.push_back(char('0' + ii));
    }
}
