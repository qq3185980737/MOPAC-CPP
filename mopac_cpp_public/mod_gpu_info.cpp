// ============================================================================
// 【中文导读】mod_gpu_info.cpp —— MOZYME 局域轨道/成键
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 mod_gpu_info.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// mod_gpu_info.cpp
#include "mod_gpu_info.h"
void getGPUInfo(int* h, int*, int* nd, char*, int*, size_t*, int*, int*,
                int*) { *h=0; *nd=0; }
void setDevice(int, int* stat) { *stat=0; }
