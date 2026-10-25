// ============================================================================
// 【中文导读】upcase.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 upcase.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// upcase.cpp
#include "upcase.h"
#include <cctype>
void upcase(std::string& s, int n) {
    std::string buf=s;
    for (int i=0;i<n;++i) {
        char c=s[i];
        if (c>='a'&&c<='z') s[i]=char(c-'a'+'A');
        else if (c=='\t') s[i]=' ';
    }
    size_t p=s.find("EXTERNAL");
    if (p!=std::string::npos) {
        size_t q=s.find(' ',p+1);
        if (q==std::string::npos) q=s.size();
        for (size_t k=p+9;k<q;++k) s[k]=buf[k];
    }
}
