// ============================================================================
// 【中文导读】ccrep.h —— 结构/成键/电荷接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// ccrep.h — C++ translation of MOPAC 2016 "ccrep.F90".
#pragma once

// ccrep: core-core repulsion between atoms ni,nj.
// Input r (Angstrom) is converted; gab is the monopole factor.
// enuclr is the output repulsion.
void ccrep(int ni, int nj, double& r, double gab, double& enuclr);
