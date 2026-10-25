// ============================================================================
// 【中文导读】timout.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 timout.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// timout.cpp
#include "timout.h"
#include "chanel_C.h"
#include "molkst_C.h"
#include <cstdio>
#include <cmath>
#include <string>
using namespace molkst_C;
extern double wall_clock_0, wall_clock_1;
extern double CPU_0, CPU_1;
static const char* plural(int n, const char* s, const char* p) { return n > 1 ? p : s; }
void timout(int nout) {
    (void)nout;
    if (keywrd.find("LOCATE-TS") != std::string::npos) return;
    double wdays = (wall_clock_1 - wall_clock_0) / 86400.0;
    int wd = (int)std::floor(wdays);
    double whours = (wdays - wd) * 24.0;
    int wh = (int)std::floor(whours);
    double wmins = (whours - wh) * 60.0;
    int wm = (int)std::floor(wmins);
    double ws = (wmins - wm) * 60.0;
    if (wd > 0) std::fprintf(chanel_C::out_unit(), "          WALL-CLOCK TIME         =%d%s %d%s %d%s AND%7.3f SECONDS\n", wd, plural(wd," DAY "," DAYS"), wh, plural(wh," HOUR "," HOURS"), wm, plural(wm," MINUTE "," MINUTES"), ws);
    else if (wh > 0) std::fprintf(chanel_C::out_unit(), "          WALL-CLOCK TIME         =%d%s %d%s AND%7.3f SECONDS\n", wh, plural(wh," HOUR "," HOURS"), wm, plural(wm," MINUTE "," MINUTES"), ws);
    else if (wm > 0) std::fprintf(chanel_C::out_unit(), "          WALL-CLOCK TIME         =%d%s AND%7.3f SECONDS\n", wm, plural(wm," MINUTE "," MINUTES"), ws);
    else std::fprintf(chanel_C::out_unit(), "          WALL-CLOCK TIME         =%15.3f SECONDS\n", ws);
    double ctim = CPU_1 - CPU_0;
    double cdays = ctim / 86400.0;
    int cd = (int)std::floor(cdays);
    double chours = (cdays - cd) * 24.0;
    int ch = (int)std::floor(chours);
    double cmins = (chours - ch) * 60.0;
    int cm = (int)std::floor(cmins);
    double cs = (cmins - cm) * 60.0;
    if (cd > 0) std::fprintf(chanel_C::out_unit(), "          COMPUTATION TIME        =%d%s %d%s %d%s AND%7.3f SECONDS\n", cd, plural(cd," DAY "," DAYS"), ch, plural(ch," HOUR "," HOURS"), cm, plural(cm," MINUTE "," MINUTES"), cs);
    else if (ch > 0) std::fprintf(chanel_C::out_unit(), "          COMPUTATION TIME        =%d%s %d%s AND%7.3f SECONDS\n", ch, plural(ch," HOUR "," HOURS"), cm, plural(cm," MINUTE "," MINUTES"), cs);
    else if (cm > 0) std::fprintf(chanel_C::out_unit(), "          COMPUTATION TIME        =%d%s AND%7.3f SECONDS\n", cm, plural(cm," MINUTE "," MINUTES"), cs);
    else std::fprintf(chanel_C::out_unit(), "          COMPUTATION TIME        =%15.3f SECONDS\n", cs);
}
