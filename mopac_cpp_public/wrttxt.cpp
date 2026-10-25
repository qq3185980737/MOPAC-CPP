// wrttxt.cpp
// ============================================================================
// 【中文导读】wrttxt() —— 输出标题/关键字回显文本
// ----------------------------------------------------------------------------
// 对应 wrttxt.F90。把 refkey（关键字回显）与 koment/title 打印到主输出；
// 若关键字里含 CHAINS/START_RES 且 refkey 尚未记录，则补到 refkey[0]。
// " NULL" 是 MOPAC 的哨兵，表示该行未设置，需要跳过。
// ============================================================================

#include "wrttxt.h"
#include "chanel_C.h"
#include <string>
#include <cstdio>
#include <cctype>
#include "molkst_C.h"
using namespace molkst_C;

// wrttxt(): 打印关键字回显/注释/标题。

void wrttxt(int iprt) {
    bool l_chains = false, l_start = false;
    for (int i = 0; i < 6; ++i) {
        if (refkey[i].find(" NULL") != std::string::npos) break;
        line = " " + refkey[i];
        for (auto& ch : line) ch = (char)toupper(ch);
        if (!l_chains) l_chains = (line.find(" CHAINS") != std::string::npos);
        if (!l_start)  l_start  = (line.find(" START_RES") != std::string::npos);
    }
    size_t i = keywrd.find(" CHAINS");
    if (i != std::string::npos && !l_chains) {
        size_t j = keywrd.find(")", i + 7) + i + 7;
        refkey[0] = keywrd.substr(i, j - i) + refkey[0];
    }
    i = keywrd.find(" START_RES");
    if (i != std::string::npos && !l_start) {
        size_t j = keywrd.find(")", i + 10) + i + 10;
        refkey[0] = keywrd.substr(i, j - i) + refkey[0];
    }
    const bool to_file = (chanel_C::out_fp != nullptr);
    auto rstrip = [](std::string s) {
        while (!s.empty() && (s.back() == ' ' || s.back() == '\t')) s.pop_back();
        return s;
    };
    for (int k = 0; k < 6; ++k) {
        if (refkey[k].find(" NULL") != std::string::npos) break;
        std::string s = to_file ? rstrip(refkey[k]) : refkey[k];
        if (to_file && (s.empty() || s == " ")) continue;   // 档案文件不写空行
        std::fprintf(chanel_C::out_unit(), "%s\n", s.c_str());
    }
    // F90 readmo writes the comment/title records with a leading blank.
    if (koment.find(" NULL ") == std::string::npos) {
        std::string s = to_file ? rstrip(koment) : koment;
        if (!(to_file && s.empty()))
            std::fprintf(chanel_C::out_unit(), " %s\n", s.c_str());
    }
    if (title.find(" NULL ") == std::string::npos) {
        std::string s = to_file ? rstrip(title) : title;
        if (!(to_file && s.empty()))
            std::fprintf(chanel_C::out_unit(), " %s\n", s.c_str());
    }
}
