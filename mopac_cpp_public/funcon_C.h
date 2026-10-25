// ============================================================================
// 【中文导读】funcon_C.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// funcon_C.h — C++ mapping of Fortran module "funcon_C".
#pragma once
namespace funcon_C {
extern double pi;
extern double a0;      // Bohr radius in Angstrom
extern double ev;       // 1 Hartree in eV
extern double fpc_9;    // 1 eV in kcal/mol
extern double fpc_10;   // Avogadro's number
extern double fpc_8;    // speed of light (cm/s)
extern double fpc_6;    // Planck constant (erg-sec)
extern double fpc_2;
extern double fpc_1;
extern double fpc_5;    // R (cal/mol/K)
extern double fpc_7;    // Boltzmann (erg/K)
}  // namespace funcon_C