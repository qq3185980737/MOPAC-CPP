// ============================================================================
// 【中文导读】genun.h —— 数值/几何工具接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// genun.h — C++ translation of MOPAC 2016 "genun.F90".
#pragma once
#include <vector>
void genun(std::vector<std::vector<double>>& u, int& n);
bool collid(double rw, const std::vector<double>& cw,
            const std::vector<std::vector<double>>& cnbr,
            const std::vector<double>& rnbr, int nnbr, int ishape);
