// ============================================================================
// 【中文导读】write_pressure.h —— 输出接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明 write_pressure(iprt)：对应 moldat.F90:1174-1248 的子程序，输出
// “约束平移向量所需的压力”（周期体系 id==3，把梯度投影到 Tv 方向再换算
// 成 GPa 写入 press(1:3)），并给出提示/警告。由 writmo 周期块调用。
// ============================================================================

// write_pressure.h — C++ translation of moldat.F90:1174-1248
// ("write_pressure").
#pragma once

// write_pressure(iprt):
//   iprt < 0  -> no output
//   iprt == 0 -> to_screen()
//   iprt > 0  -> printf() to stdout (list-directed style adds one leading blank)
void write_pressure(int iprt);
