// to_screen.cpp — screen/log messages (F90: write(0,...)); keep them off the
// .out file, which is stdout.

//
// ============================================================================
// 【中文导读】to_screen() —— 屏幕/日志消息输出
// ----------------------------------------------------------------------------
// 对应 F90 的 write(0,...)。iw0 是 GUI/屏幕摘要通道：iw0 > -1 时消息写到
// stderr（不进入 .out），这样控制台能看到进度，而 .out 文件保持纯结果。
// ============================================================================

#include "to_screen.h"
#include "chanel_C.h"
#include <cstdio>

// to_screen(): 仅当 iw0 > -1 时把消息写到 stderr。

void to_screen(const std::string& text) {
    if (chanel_C::iw0 > -1) { std::fprintf(stderr, "%s\n", text.c_str()); }
}
