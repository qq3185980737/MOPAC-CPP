// ============================================================================
// 【中文导读】ligand.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// ligand.h — C++ translation of "ligand.F90".
#pragma once
#include <string>
#include <vector>

// ligand: re-labels atoms/residues in a PDB-like input: detects phosphate,
// sulfate, water, ethylene glycol, glycerol, and general hetero groups; applies
// the " XENO(...)" keyword overrides; writes "ATOM  "/"HETATM" records into
// txtatm.  ires/nfrag advance as residues are consumed.
void ligand(int& ires, const std::vector<int>& start_res, int& nfrag);

// moiety: collects all atoms connected (through bonds) to atom istart into the
// "used" list; hydrogen atoms are moved to the end of the list.  On return
// lused(new+1..) holds the moiety atoms.
void moiety(std::vector<bool>& iopt, std::vector<int>& lused, int istart,
            int& n_new);

// nheavy: number of non-hydrogen atoms bonded to atom icc.
int nheavy(int icc);

// identify_hexose: determine which hexose is present from the chiral centers.
void identify_hexose(int ninres, const std::vector<int>& inres,
                     std::string& nam, std::string& name);

// inc_res: advance the residue counter, honoring the start_res fragment list.
void inc_res(int& ires, const std::vector<int>& start_res, int& nfrag);
