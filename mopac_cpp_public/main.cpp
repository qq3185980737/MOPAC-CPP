// main.cpp — C++ translation of "MOPAC for Windows.F90" (program MOPAC_win).
// The original uses Intel QuickWin (IFQWIN/DFLIB) for a windowed menu GUI.
// That platform glue is dropped; this provides a console entry point that
// performs the same core initialization and runs MOPAC.

//
// ============================================================================
// 【中文导读】main() —— 控制台入口（原版为 QuickWin 窗口程序 MOPAC_win）
// ----------------------------------------------------------------------------
// 1) 调用链：main() -> run_mopac() -> getdat() 读取 <job>.mop -> 按关键字分派
//    计算（iter/ef/lbfgs/drc/...）-> writmo() 输出结果。
// 2) 输入文件来自命令行：main() 把 run 置为 1；getdat() 用 __argv[run] 即
//    argv[1] 作为数据卡文件名，并据此拼出 .out/.arc 等文件名。
// 3) gui=false 表示控制台模式；iw0 是 GUI 用的“摘要输出”通道（原版默认 -1，
//    表示不用），这里置 0。真正的输出走 stdout，会被 run_mopac() 重定向到
//    <input>.out。
// 4) 许可证：原版通过 password() 校验并设置 ijulian/site_no；本归档的
//    mopac.F90 中该调用已被注释掉。本移植直接赋固定值（ijulian=123,
//    site_no=88），等价于跳过校验，
//    仅供研究与本地验证，不要用于对外分发。
// 5) 注意：argv 形参本身没有被使用，取文件名发生在 getdat() 内部（见上）。
// ============================================================================

#include <cstdio>
#include <cstdlib>
#include <string>

#include "chanel_C.h"
#include "molkst_C.h"

// Forward declarations of the computational entry points.
void GetDateStamp(std::string& line, std::string& verson);
void run_mopac();

using namespace chanel_C;
using namespace molkst_C;

int main(int argc, char** argv) {
    gui = false;
    iw0 = 0;

    // 许可证字段：原版由 password() 设置；这里给固定值以跳过校验（见文件头说明）。

    ijulian = 123;
    site_no = 999;              // match reference build banner ("Site#:  999")
    run = 1;                    // argv[1] = data-set filename
    program_name = "Standalone MOPAC ";

    if (argc <= 1) {
        std::printf("\n");
        std::printf("  MOPAC: supply a data-set filename as the first argument.\n");
        return EXIT_SUCCESS;
    }

    GetDateStamp(line, verson);
    run_mopac();
    return EXIT_SUCCESS;
}
