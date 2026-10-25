// ============================================================================
// 【中文导读】fdate.cpp —— 主驱动/运行控制
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 fdate.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// fdate.cpp — Fortran fdate intrinsic: current date/time string (writmo dep).
#include "fdate.h"
#include <ctime>
#include <string>

void fdate(std::string& date_str) {
  std::time_t t = std::time(nullptr);
  std::tm tm{};
  localtime_s(&tm, &t);
  char buf[40];
  std::strftime(buf, sizeof(buf), "%a %b %e %H:%M:%S %Y", &tm);
  date_str = buf;
}
