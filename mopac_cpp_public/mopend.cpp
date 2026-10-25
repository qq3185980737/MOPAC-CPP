// ============================================================================
// 【中文导读】mopend.cpp —— MOZYME 局域轨道/成键
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mopend.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mopend.cpp — C++ translation of mopend.F90 (mopend + summary).
#include "mopend.h"

#include <cstdio>
#include <string>
#include <vector>

#include "molkst_C.h"

using namespace molkst_C;

void summary(const std::string& txt, int ntxt);

void mopend(const char* txt) { mopend(std::string(txt)); }
void mopend(const std::string& txt) {
    moperr = true;
    errtxt = txt;
    summary(txt, (int)txt.size());
    if (txt != "JOB ENDED NORMALLY")
        std::printf("\n%10s%s\n", "", txt.c_str());
}

void summary(const std::string& txt, int ntxt) {
    static int nmessages = 0;
    static bool first = true;
    static std::vector<std::string> messages(20, "");
    const int lim = 120;
    if (first) {
        messages[0] = "JOB ENDED NORMALLY";
        first = false;
    }
    if (ntxt == 1) {
        // Print the termination-message box.
        int max_txt = 1;
        for (int i = 0; i < nmessages; ++i)
            max_txt = std::max(max_txt, (int)messages[i].size());
        max_txt = std::min(lim, std::max(max_txt + 4, 22));
        bool bad = messages[0].substr(0, 18) != "JOB ENDED NORMALLY";
        if (bad) max_txt = std::max(max_txt, 78);

        std::string top(max_txt, '*');
        std::string blank(max_txt, ' ');
        blank[0] = '*';
        blank[max_txt - 1] = '*';

        std::printf("\n %s\n", top.c_str());
        std::printf(" %s\n", blank.c_str());
        if (bad) {
            std::string line = "*     Error and normal termination messages reported in this calculation";
            line += std::string(max_txt - (int)line.size(), ' ');
            line[max_txt - 1] = '*';
            std::printf(" %s\n", line.c_str());
            std::printf(" %s\n", blank.c_str());
            for (int i = 0; i < nmessages; ++i) {
                if (messages[i].find("JOB ENDED NORMALLY") != std::string::npos) continue;
                line = "* " + messages[i];
                if ((int)line.size() < max_txt)
                    line += std::string(max_txt - (int)line.size(), ' ');
                line[max_txt - 1] = '*';
                std::printf(" %s\n", line.c_str());
            }
        }
        std::string line = "* JOB ENDED NORMALLY";
        if ((int)line.size() < max_txt)
            line += std::string(max_txt - (int)line.size(), ' ');
        line[max_txt - 1] = '*';
        std::printf(" %s\n", line.c_str());
        std::printf(" %s\n", blank.c_str());
        std::printf(" %s\n", top.c_str());

        nmessages = 0;
        messages[0] = "JOB ENDED NORMALLY";
    } else {
        if (nmessages == 20) return;
        ++nmessages;
        int i = std::min(ntxt, lim - 3);
        messages[nmessages - 1] = txt.substr(0, i);
    }
}
