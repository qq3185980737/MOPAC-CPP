// ============================================================================
// 【中文导读】convert_storage.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// convert_storage.h — C++ translation of MOPAC 2016 "convert_storage.F90".
#pragma once
#include <vector>

// MOZYME packed -> lower-triangle packed matrix.
void convert_mat_packed_to_triangle(const std::vector<double>& matrix_packed,
                                    std::vector<double>& matrix_triangle);

// MOZYME packed LMOs -> conventional square eigenvector matrix.
void convert_lmo_packed_to_square(std::vector<std::vector<double>>& c_square);
