// ============================================================================
// 【中文导读】gpu_info.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// gpu_info.h — C++ binding for the external GPU info / device-selection C API.
// Mirrors mod_gpu_info.F90 (getGPUInfo / setDevice).
#pragma once
#include <cstddef>
#include <cstdint>

extern "C" {
// Query GPU availability. Arrays are sized up to 6 devices.
void getGPUInfo(bool* hasGpu, bool* hasDouble, int* nDevices,
                char* name, int* name_size, std::size_t* totalMem,
                int* clockRate, int* major, int* minor);

// Select the active CUDA device.
void setGPU(int idevice, bool* stat);
}
