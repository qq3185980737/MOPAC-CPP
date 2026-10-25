// ============================================================================
// 【中文导读】adjvec.h —— MOZYME 接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// adjvec.h — C++ translation of MOPAC 2016 "adjvec.F90".
//
// Orthogonalise one localised molecular orbital (LMO) vector against another
// in the MOZYME sparse (CSR-like) packed format.  Arrays keep their Fortran
// 1-based indexing (index 0 unused / padding).
#pragma once

#include <vector>

void adjvec(std::vector<double>& cvecb, int ncvb,
            std::vector<int>& icvecb, int nib,
            const std::vector<int>& nncb, std::vector<int>& ncb_loc, int nnb,
            const std::vector<int>& ncvecb, int lmob,
            const std::vector<int>& iorbs,
            const std::vector<double>& cveca, int ncva,
            const std::vector<int>& icveca, int nia,
            const std::vector<int>& nnca, const std::vector<int>& nca_loc,
            int nna_loc, const std::vector<int>& ncveca, int lmoa,
            double beta, std::vector<int>& iused, double& sumtot);
