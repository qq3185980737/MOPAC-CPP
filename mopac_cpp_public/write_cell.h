// ============================================================================
// 【中文导读】write_cell.h —— 输出接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// write_cell.h — C++ translation of "write_cell" (moldat.F90).
#pragma once
// write_cell(iprt): print unit-cell data and formula-unit properties.
// iprt < 0 or GUI: no output; iprt == 0: only print if nstep changed.
void write_cell(int iprt);
