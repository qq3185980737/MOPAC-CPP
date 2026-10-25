// ============================================================================
// 【中文导读】mod_vars_cuda.h —— MOZYME 接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// mod_vars_cuda.h — C++ translation.
#pragma once
namespace mod_vars_cuda {
    constexpr int nthreads_gpu = 256;
    constexpr int nblocks_gpu = 256;
    extern bool lgpu;
    constexpr int prec = 8;
    extern int ngpus, gpu_id;
}
