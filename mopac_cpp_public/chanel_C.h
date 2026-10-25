// chanel_C.h — C++ mapping of Fortran module "chanel_C".
//
// ============================================================================
// 【中文导读】chanel_C —— 文件通道号与文件名（对应 F90 module chanel_C）
// ----------------------------------------------------------------------------
// 原版用 iw/ir/iarc/... 作为 Fortran unit 号，围绕它们 open/read/write/close。
// C++ 版保留了这些符号，但语义已经变化：
//   - 主输出：不再用 iw 写文件；run_mopac() 把 stdout 重定向到 <job>.out，
//     子程序用 printf/fprintf(stdout) 输出（iw 值只作占位/兼容参数）；
//   - 输入/档案/密度等文件：由 ifstream/ofstream 按 *_fn 文件名直接打开；
//   - open/close/rewind 在 C++ 版中多为 no-op，代码里保留原版注释作对照。
//
// 重要：C++ 侧的 unit 数值并不严格沿用 Fortran（chanel_C.cpp 中 iw=6、ir=5、
// iarc=8、ires=9、iden=4、iscr=12、ilog=7、lbfgs_it=13），不要按 Fortran
// 手册的 unit 约定推断；以实际调用点为准。
//
// 文件名变量（*_fn）在 getdat()/setup 阶段按作业名拼出，例如 output_fn、
// archive_fn、density_fn 等；新增/迁移输出时优先使用这些变量。
// ============================================================================

#pragma once
#include <cstdio>
#include <string>
namespace chanel_C {
extern int iw; extern int iw0; extern int ires; extern int iarc, ixyz;     // main output unit
// out_fp: 可选的输出重定向流。为 nullptr 时所有输出走 stdout（默认）；
// 写 archive/summary 文件时由 writmo 临时指向该文件，让 geout/wrttxt/
// dimens/timout 等既有的 stdout 输出函数复用同一套格式写到文件。
extern std::FILE* out_fp;
inline std::FILE* out_unit() { return out_fp != nullptr ? out_fp : stdout; }
extern int iump;  // grid surface data unit (grid.F90)
extern int ir;
extern int imep;    // input data unit (getsym.F90)
extern int iden;  // density file unit (pinout.F90)
extern int igpt;   // GRAPH file unit (mullik.F90)
extern bool log;   // true if a log file is open
extern int ilog;   // log output unit
extern std::string input_fn, output_fn, restart_fn, density_fn, log_fn, end_fn;
extern int lbfgs_it;  // iteration-summary unit (mainlb)
extern std::string archive_fn, brillouin_fn, esp_fn, ump_fn, mep_fn, pol_fn;
extern std::string gpt_fn, esr_fn, xyz_fn, syb_fn, cosmo_fn;
extern std::string job_fn;   // current job / dataset filename
extern int iscr;   // scratch file unit (geoutg.F90)
}
