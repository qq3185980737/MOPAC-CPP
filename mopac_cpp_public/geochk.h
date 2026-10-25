// ============================================================================
// 【中文导读】geochk.h —— 输入/几何接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// geochk.h — C++ translation of MOPAC 2016 "geochk.F90" (3493 lines,
// 12 subprograms: geochk, extvdw_for_MOZYME, fix_charges, add_sp_H,
// add_sp2_H, add_sp3_H, compare_sequence, update_txtatm, rectify_sequence,
// write_sequence, site, find_salt_bridges).
#pragma once
#include <string>
#include <vector>

// Main driver: checks the geometry / Lewis structure, identifies ionized
// atoms, computes the system charge, optionally resequences residues.
void geochk();

// extvdw_for_MOZYME lives in extvdw_for_MOZYME.cpp.
// find_salt_bridges lives in find_salt_bridges.cpp.

// Modify CHARGE=n keywords in refkey(1) and keywrd so the system runs with MOZYME.
void fix_charges(int ichrge);

// Place a hydrogen atom at the apex of a triangle given by i1-i-i2 (sp).
void add_sp_H(int i1, int i, int i2);
// Place a hydrogen atom in the plane of i1-i-i2 (sp2).
void add_sp2_H(int i1, int i, int i2);
// Place a hydrogen atom at the apex of a tetrahedron i1-i-i2-i3 (sp3).
void add_sp3_H(int i1, int i, int i2, int i3);

// Compare calculated residue names with those in the original data-set
// (txtatm1) and print any differences.
void compare_sequence(int n_new);

// Re-label atoms in txtatm (hydrogen numbering, chain letters, PDB numbers).
void update_txtatm(bool output, bool sort);

// Move out-of-sequence atoms to their correct position after a chain break.
void rectify_sequence();

// Write out all residue names (one-letter codes) using txtatm.
void write_sequence();

// Ionize / de-ionize residues and add or remove hydrogen atoms as requested
// by keyword SITE.
void site(const std::vector<bool>& neutral, const std::vector<char>& chain,
          const std::vector<int>& res, std::vector<std::vector<char>>& charge,
          int nres, int max_sites, std::string& allkey);
