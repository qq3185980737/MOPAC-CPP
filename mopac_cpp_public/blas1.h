// ============================================================================
// 【中文导读】blas1.h —— 数值/几何工具接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// blas1.h — standard BLAS Level-1 routines used by MOPAC.
#pragma once
#include <vector>

// y = a*x + y
void daxpy(int n, double a, const std::vector<double>& x, int incx,
           std::vector<double>& y, int incy);
// x = a*x
void dscal(int n, double a, std::vector<double>& x, int incx);
// swap x and y
void dswap(int n, std::vector<double>& x, int incx,
           std::vector<double>& y, int incy);
// dot product
double ddot(int n, const std::vector<double>& x, int incx,
            const std::vector<double>& y, int incy);
// Euclidean norm
double dnrm2(int n, const std::vector<double>& x, int incx);
// copy x -> y
void dcopy(int n, const std::vector<double>& x, int incx,
           std::vector<double>& y, int incy);
// index of largest absolute value (1-based result)
int idamax(int n, const std::vector<double>& x, int incx);
