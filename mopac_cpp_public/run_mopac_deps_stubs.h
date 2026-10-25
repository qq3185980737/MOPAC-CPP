// ============================================================================
// 【中文导读】run_mopac_deps_stubs.h —— 主驱动/运行控制接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// run_mopac_deps_stubs.h — declarations for routines called by run_mopac.
// All are implemented in their real translation files (see .cpp); this header
// exists so old-style callers can include one header for the whole set.
#pragma once
#include <string>

// update_txtatm: replace input atom labels with PDB-style labels (geochk.F90,
// translated in geochk_txtatm.cpp).
void update_txtatm(bool output, bool sort);

// write_sequence: print the residue sequence (geochk.F90, geochk_txtatm.cpp).
void write_sequence();

// delete_MOZYME_arrays: free MOZYME local-orbital arrays
// (set_up_MOZYME_arrays.F90, set_up_MOZYME_arrays.cpp).
void delete_MOZYME_arrays();

// Locate_TS_for_Proteins: locate a transition state for proteins (big_swap.cpp).
void Locate_TS_for_Proteins();

// Refine_TS_for_Proteins: refine a transition state for proteins (big_swap.cpp).
void Refine_TS_for_Proteins();

// MKL thread control (real MOPAC calls MKL; minimal stand-in here).
int mkl_get_max_threads();
void mkl_set_num_threads(int n);

// write_path_html: HTML output of a reaction path (pathk.F90, pathk.cpp).
void write_path_html();

// add_path: prepend the current working directory to a relative filename
// (readmo.F90, translated in pathk.cpp / pdbout.cpp).
void add_path(std::string& fn);
