// ============================================================================
// 【中文导读】myword.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 myword.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// myword.cpp — C++ translation of MOPAC 2016 "myword.F90".
// F90: j = index(keywrd,testwd); do while keywrd(j:j)==' ': j=j+1  (check-then-advance).
#include "myword.h"
bool myword(std::string& keywrd, const std::string& testwd) {
    bool res=false;
    size_t len_key=keywrd.find_last_not_of(' ');
    if (len_key==std::string::npos) return false;
    bool quote=false;
    size_t j;
    while ((j=keywrd.find(testwd))!=std::string::npos) {
        while (j<keywrd.size() && keywrd[j]==' ') j++;
        res=true;
        size_t k;
        for (k=j;k<=len_key;++k) {
            if (keywrd[k]=='"') quote=!quote;
            if (!quote) {
                if (keywrd[k]=='=' || keywrd[k]==' ') {
                    if (keywrd[k]=='=') { keywrd[k]=' '; continue; }
                    size_t j2;
                    for (j2=k+1;j2<=len_key;++j2) {
                        if (keywrd[j2]=='=') { keywrd[j2]=' '; break; }
                        if (keywrd[j2]!=' ') break;
                    }
                    break;
                }
            }
            keywrd[k]=' ';
        }
    }
    return res;
}
