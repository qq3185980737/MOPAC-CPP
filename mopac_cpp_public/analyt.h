// ============================================================================
// 【中文导读】analyt.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// analyt.h — C++ translation of MOPAC 2016 "analyt.F90".
#pragma once
#include <vector>

// analyt: analytical derivatives between two atoms i=2, j=1.
//   psum/palpha/pbeta: packed 2D density / Pi arrays (size mpack+1).
//   coord[dim][atom] (atom 1..2). nat[1..2] = atomic numbers.
//   eng[1..3] output derivative w.r.t. each Cartesian component.
void analyt(const std::vector<double>& psum, const std::vector<double>& palpha,
            const std::vector<double>& pbeta,
            const std::vector<std::vector<double>>& coord,
            const int nat[3], int jja, int jjd, int iia, int iid,
            double eng[4]);

// Internal (also called from Fortran as separate routines).
void delmol(const std::vector<std::vector<double>>& coord, int i, int j,
            int ni, int nj, int ia, int id, int ja, int jd, int ix,
            double rij, double tomb, int& isp);
void delri(double dg[23], int ni, int nj, double rr, double del1);
void rotat(const std::vector<std::vector<double>>& coord, int i, int j, int ix,
           int idx, double rij, double del1);
void ders(int m, int n, double rr, double del1, double del2, double del3,
          int is, int iol);
