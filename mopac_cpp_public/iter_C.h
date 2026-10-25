// iter_C.h — C++ mapping of Fortran module "iter_C" (subset used).
//
// ============================================================================
// 【中文导读】iter_C —— SCF 迭代工作数组（对应 F90 module iter_C）
// ----------------------------------------------------------------------------
// 这些数组只被 iter()/pulay()/interp() 等迭代相关代码使用，保存“上一轮/上
// 两轮”的密度和收敛加速器的中间结果：
//   - pold/pold2/pold3、pbold/pbold2/pbold3：密度历史值（下三角打包，1-based），
//     供 cnvg() 计算变化量、pulay() 做 DIIS/Camp-King 外推；
//   - vec_ai/fock_ai/p_ai/h_ai/vecl_ai（及 _bi 版本）：Camp-King 收敛器的
//     scratch 区，大小约 (norbs+1)^2，按平铺数组使用；
//   - eold_alpha/eold_beta、theta：上一轮能量与占据数工作向量；
//   - md_test/is_PARAM/lgpu：模式标志。md_test 表示 SCF 不收敛时是否尝试
//     Camp-King 重试；is_PARAM 是 PARAM 参数化程序模式（正常 MOPAC 为 false）。
// 注意：0 号元素是 padding；这些是跨迭代/跨作业保持的全局量，不要当作局部变量。
// ============================================================================

#pragma once
#include <vector>
namespace iter_C {
// pold(:): previous-iteration density/packed vector (1-based padding).
extern std::vector<double> pold;
extern std::vector<double> pold2, pold3;
extern std::vector<double> pbold, pbold2, pbold3;
// Camp-King converger scratch (norbs x norbs, norbs*norbs)
extern std::vector<double> vec_ai, fock_ai, p_ai, h_ai, vecl_ai;
extern std::vector<double> vec_bi, fock_bi, p_bi, h_bi, vecl_bi;
extern double eold_alpha, eold_beta;
extern std::vector<double> theta;
extern bool md_test, is_PARAM;
extern bool lgpu;}
