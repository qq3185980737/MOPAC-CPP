// ============================================================================
// 【中文导读】md_driver.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// md_driver.h — C++ port of the user's Fortran molecular-dynamics driver
// (originally "run_md_velocity_verlet" + "gen_md_vel" from md.txt).
//
// Runs a velocity-Verlet MD for the current Cartesian system (H2/O2 mixtures
// and beyond) using MOPAC PM7 gradients (kcal/mol/A) as forces. All inputs are
// read from the existing globals (xparam, numat, nvar, ams, loc, nat, escf,
// elect/enuclr/solv_energy); temperature / box / timestep are controlled by
// keywords TEMP=, BOX= and DT=.
#pragma once
#include <string>

void run_md_velocity_verlet(const std::string& keywrd, int iw);
