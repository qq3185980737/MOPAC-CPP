// ============================================================================
// 【中文导读】nuchar.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 nuchar.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// nuchar.cpp — C++ translation of MOPAC 2016 "nuchar.F90".
// Determines and returns the real values of all numbers found in 'line'.
#include "nuchar.h"
#include "reada.h"

#include <string>

void nuchar(char* line, int l_line, double* value, int& nvalue) {
    int istart[41];
    int i;
    bool leadsp;
    char tab, comma, space;
    tab = char(9);
    comma = ',';
    space = ' ';
    // CLEAN OUT TABS AND COMMAS
    for (i = 1; i <= l_line; ++i) {
        if (line[i - 1] != tab && line[i - 1] != comma) continue;
        line[i - 1] = space;
    }
    // FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES
    // FOLLOWED BY A CHARACTER
    leadsp = true;
    nvalue = 0;
    for (i = 1; i <= l_line; ++i) {
        if (leadsp && line[i - 1] != space) {
            ++nvalue;
            istart[nvalue] = i;
        }
        leadsp = (line[i - 1] == space);
    }
    // FILL NUMBER ARRAY
    for (i = 1; i <= nvalue; ++i)
        value[i] = reada(std::string(line), istart[i]);
}
