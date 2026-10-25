// ============================================================================
// 【中文导读】write_unit_cell_hof.h —— 输出接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明 write_unit_cell_HOF(iprt)：对应 moldat.F90:1118-1173 的子程序，
// 打印“每个单胞的生成热（H.o.F. per unit cell）”及参与计算的单胞数。
// 仅在周期体系（id==3）且 writmo 主输出路径中调用。
// ============================================================================

// write_unit_cell_hof.h — C++ translation of moldat.F90
// "write_unit_cell_HOF" (lines 1118-1173).
#pragma once

// write_unit_cell_HOF(iprt):
//   iprt < 0        -> no output (disabled stream)
//   iprt == 0       -> to_screen() only
//   iprt > 0        -> printf() to stdout (MOPAC's iw is stdout)
// iprt 为 0 时不会写 .out 文件；writmo 主路径传 iw=6。
void write_unit_cell_HOF(int iprt);
