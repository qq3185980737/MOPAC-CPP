// ============================================================================
// 【中文导读】eigen.h —— 对角化/本征值接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// eigen.h — C++ translation of MOPAC 2016 "eigen.F90".
#pragma once
#include <vector>

// Parse VECTORS keyword -> print ranges.
void eigen_limits(int& print_nocc, int& print_nvir);
void eigen(bool write_gpt, bool write_out);

// MOZYME: build the Fock matrix in the LMO basis (eigen.F90:192-379).
void FlmoFromFao(int nlmo, const std::vector<int>& ncx, const std::vector<int>& nncx,
                 const std::vector<int>& ncxxx, int n_dim,
                 std::vector<char>& latoms, std::vector<double>& ws,
                 const std::vector<int>& icxxx, int icxxx_dim,
                 const std::vector<double>& cxxx, int cxxx_dim,
                 std::vector<double>& ac, std::vector<double>& fmo_expanded);

// MOZYME: diagonalise the LMO-basis Fock matrix and expand back to AO basis
// (eigen.F90:438-499). Offsets support the virtual set starting at noccupied+1.
void lmo_to_eigenvectors(int nlmo, const std::vector<int>& ncx, const std::vector<int>& nncx,
                         const std::vector<int>& ncxxx, int n_dim,
                         const std::vector<int>& icxxx, int icxxx_dim,
                         const std::vector<double>& cxxx, int cxxx_dim,
                         std::vector<double>& eigs_out, int eig_offset,
                         std::vector<std::vector<double>>& c, int col_offset);
