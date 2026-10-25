// ============================================================================
// 【中文导读】analyt_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// analyt_C.h — C++ mapping of Fortran module "analyt_C".
#pragma once
namespace analyt_C {
constexpr int NZTYPE_MAX = 107;
extern int nztype[NZTYPE_MAX + 1];   // 1-based
extern double ds[17];                // dimension(16)
extern double dg[23];                // dimension(22)
extern double dr[101];               // dimension(100)
extern double g[23];                 // dimension(22)
extern double tx[4], ty[4], tz[4];         // dimension(3)
extern double tdx[4], tdy[4], tdz[4];      // dimension(3)
}
