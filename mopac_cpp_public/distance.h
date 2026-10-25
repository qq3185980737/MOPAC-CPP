// ============================================================================
// 【中文导读】distance.h —— 结构/成键/电荷接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// distance.h — C++ translation of distance/angle/torsion functions from
// H_bond_correction_bits.F90 (lines 267-299).
#pragma once
#include <vector>

// Interatomic distance between atoms a and b, with lattice translation
// (id != 0) taken into account.
double distance(int a, int b);

// Bond angle at b between a and c (radians).
double angle(int a, int b, int c);

// Dihedral angle between planes (i,j,k) and (j,k,l) (radians).
double torsion(int i, int j, int k, int l);
