// ============================================================================
// 【中文导读】getval.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 getval.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// getval.cpp — C++ translation of MOPAC 2016 "getval.F90".
// F90: numeric iff ch1 AND ch2 are both NOT uppercase 'A'..'Z'; x = reada(line,1).
#include "getval.h"
#include "reada.h"
#include <string>

void getval(const std::string& line, double& x, std::string& t) {
    char ch1 = line.size() > 0 ? line[0] : ' ';
    char ch2 = line.size() > 1 ? line[1] : ' ';
    bool letter1 = (ch1 >= 'A' && ch1 <= 'Z');
    bool letter2 = (ch2 >= 'A' && ch2 <= 'Z');
    if (!letter1 && !letter2) {
        x = reada(line, 1);
        t = " ";
    } else {
        size_t sp = line.find(' ');
        t = line.substr(0, sp);
        x = -999.0;
    }
}
