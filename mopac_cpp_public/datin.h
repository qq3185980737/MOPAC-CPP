// ============================================================================
// 【中文导读】datin.h —— 结构/成键/电荷接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// datin.h — C++ translation of MOPAC 2016 "datin.F90".
#pragma once
#include <string>
#include <vector>

// Parse the EXTERNAL=/PARAMS= file list out of the keyword string.
std::vector<std::string> split_param_files(const std::string& keywrd);

// Element symbol by atomic number (1..107), matching MOPAC elemnt table.
std::string elemnt_sym(int z);

// Full parameter-file read-in (file I/O; body reads files and calls update()).
void datin(int iw);
