// ============================================================================
// 【中文导读】set.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// set.h — C++ translation.
#pragma once
namespace overlaps_C {
    extern thread_local int isp, ips;
    extern thread_local double sa, sb;
    extern thread_local double a[64], b[64];
    extern double fact[18];
}
void set_fn(double s1, double s2, int na, int nb, double rab, int ii);
