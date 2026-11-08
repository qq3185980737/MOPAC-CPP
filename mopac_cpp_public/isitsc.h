// ============================================================================
// 【中文导读】isitsc.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// isitsc.h — C++ translation.
#pragma once
void isitsc(double escf, double selcon, double emin, int& iemin,
            int& iemax, bool& okscf, int niter, int itrmax);
// MOZYME+MD conservation fix: clears the file-scope SCF plateau history
// (scf1/escf0[]) that isitsc keeps between calls.  Must be called from the
// per-step reset block in iter_for_MOZYME (icalcn != step_num) so that every
// MD frame starts SCF from a clean state instead of borrowing "plateau
// credit" from a previous frame's converged energy.
void isitsc_reset();
