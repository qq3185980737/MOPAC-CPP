// ============================================================================
// 【中文导读】mod_iso_c_binding.h —— MOZYME 接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// mod_iso_c_binding.h — C++ translation.
#pragma once
#include <cstdint>
namespace iso_c_binding {
    using c_intptr_t = intptr_t;
    constexpr int c_char = 1;
    constexpr int c_bool = 1;
    constexpr int c_int = 4;
    constexpr int c_float = 4;
    constexpr int c_long = 4;
    constexpr int c_double = 8;
    using c_size_t = size_t;
    struct c_ptr { c_intptr_t ptr; };
}
