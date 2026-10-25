// ============================================================================
// 【中文导读】delsta.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// delsta.h — C++ translation of MOPAC 2016 "delsta.F90".
#pragma once
#include <vector>

// Point-charge electrostatic derivative between atoms ii,jj.
void delsta(const std::vector<int>& nat, const std::vector<int>& iorbs,
            const std::vector<double>& p,
            const std::vector<std::vector<double>>& cdi,
            std::vector<double>& dstat, int ii, int jj);
