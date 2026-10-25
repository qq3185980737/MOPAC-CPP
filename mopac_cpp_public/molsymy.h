// ============================================================================
// 【中文导读】molsymy.h —— MOZYME 接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// molsymy.h — C++ translation.
#pragma once
void molsym(double* coord, int& ierror, double* r);
void chi(double toler, double* coord, int ioper, int& iqual);
void makopr(int numat, double* coord, int& ierror, double* r);
void orient(int numat, double* coord, double* r);
void plato(double* coord, double* r, bool& cubic);
bool symdec(int n1, const int* ielem);
void cartab();
