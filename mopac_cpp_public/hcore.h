// ============================================================================
// 【中文导读】hcore.h —— SCF/Fock 接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// hcore.h — C++ translation of MOPAC 2016 "hcore.F90".
#pragma once
void hcore();
// wstore: fill the one-centre two-electron integral block for atom `ni` into
// the packed W array at offset `kr` (column-major ilim x ilim), then advance
// `kr` by ilim*ilim.  For d orbitals the 243-entry repd table is used.
// 中文：填充原子 ni 的一中心双电子积分块（s/p/d），按列主序写入 W[kr...]，
// 最后 kr += ilim*ilim；不按 cutof2 丢弃（旧注释描述的行为并不存在）。
void wstore(double* w, int& kr, int ni, int ilim);
