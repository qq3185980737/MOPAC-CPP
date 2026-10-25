// ============================================================================
// 【中文导读】GetDateStamp.h —— 输入/几何接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// GetDateStamp.h — C++ translation of MOPAC 2016 "GetDateStamp.F90".
#pragma once
#include <string>
void GetDateStamp(std::string& date_time, std::string& version);
