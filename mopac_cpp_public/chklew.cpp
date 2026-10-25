// chklew.cpp — C++ translation of MOPAC 2016 "chklew.F90".
//
// ============================================================================
// 【中文导读】chklew() —— Lewis 结构统计与错误上报
// ----------------------------------------------------------------------------
// 对应 chklew.F90。维护 Lewis 电子/键计数表（add_Lewis_element 等），并在
// 电子数/键级不一致时调用 mopend/memory_error；2026-10 已删除本文件里的
// 空 mopend/memory_error stub，恢复真实错误传播。
//
// 移植进度（2026-10-25，逐段编译 + 回归）：
//   F90 128-618  主流程：σ 骨架、孤对、五元环 π 键、外层 DO 循环、开链 π 键、
//                Rule 1/2/3 重扫、环 π 键（1030/1040/1060 标签）——按原版
//                控制流逐行重译（第二轮，替换早先的"展平"写法）
//   F90 620-643  chklew 尾部：芳香碳计数分支、DEBUG/LARGE 数组转储、L 赋值
//   F90 644-652  标签 1060 死代码（仅由 F90:379 的提前出口到达，L 恒为 0）
//   F90 688-748  ring5()      五元环标记
//   F90 749-803  arom()       芳环第三个 pi 键判别
//   F90 804-899  arom2()      芳环第二个 pi 键判别
//   F90 900-975  add_Lewis_element()  Lewis 元素登记与电子/轨道记账
// 上述实现均与 Fortran 逐行对照；控制流用 C++ 的显式状态/提前出口等价表达，
// 注释保留原始 F90 行号与标签，便于审计。
// ============================================================================

#include "chklew.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"
#include "mopend.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;
using namespace MOZYME_C;
using namespace parameters_C;

// memory_error 的真实实现（std::string 版）在 setup_mopac_arrays.cpp。

namespace {

// Fortran list-directed WRITE with an "(20I4)" format: the values are printed
// right-aligned in 4-character fields, 20 values per record, and the format is
// reverted (a new record started) after every 20 values.  `pick` yields the
// value for atom number i.
template <typename Pick>
void print_i4_row(int lo, int hi, Pick pick) {
    int col = 0;
    for (int i = lo; i <= hi; ++i) {
        std::printf("%4d", pick(i));
        if (++col == 20) {
            std::printf("\n");
            col = 0;
        }
    }
    if (col != 0) std::printf("\n");
}

// F90 `write (iw, *) "literal"` -list-directed output of a character constant
// starts one column in from the left margin (the record begins with a blank
// separator), so a literal whose own first character is a blank appears with
// two leading blanks in the output file.  Reproduce that here.
void print_ld_text(const char* literal) { std::printf(" %s\n", literal); }

}  // namespace

// F90 900-975 -add_Lewis_element
// Record one Lewis structural element (bond / lone pair / virtual lone pair)
// and update the electron (IZ) and atomic-orbital (IB) bookkeeping.
void add_Lewis_element(int atom_i, int atom_j, int charge, int& element_type) {
    Lewis_tot++;
    Lewis_elem[1][Lewis_tot] = atom_i;
    Lewis_elem[2][Lewis_tot] = atom_j;
    if (atom_i > 0 && atom_j > 0) {
        iz[atom_i]--; iz[atom_j]--;
        ib[atom_i]--; ib[atom_j]--;
        element_type++;
    } else if (atom_i != 0) {
        if (atom_i > 0) {
            if (charge == -1) iz[atom_i]--;
            else if (charge == 0) iz[atom_i] -= 2;
            ib[atom_i]--;
            element_type++;
        } else {
            ib[-atom_i]--;
        }
    } else {
        if (charge == 2) iz[atom_j] -= 2;
        else if (charge == 1) iz[atom_j]--;
        ib[atom_j]--;
    }
    if (charge == 0) return;
    ions[atom_i + atom_j] += charge;
}

// F90 688-748 -ring5
// Flag the atoms of each five-membered ring through IR5 (all five atoms of a
// ring share the same non-zero flag).  The odd-looking `mb[j] < 3` test in the
// K loop is intentional: the Fortran original re-tests MB(J) there (instead of
// MB(K)), and the translation preserves that behaviour.
void ring5(int i, const std::vector<int>& mb, std::vector<int>& ir5) {
    int ni = nbonds[i];
    int nrings = 1;
    for (int j1 = 1; j1 <= ni; ++j1) {
        int j = ibonds[j1][i];
        if (mb[j] < 3) continue;
        int nj = nbonds[j];
        for (int k1 = j1 + 1; k1 <= ni; ++k1) {
            int k = ibonds[k1][i];
            if (mb[j] < 3) continue;
            int nk = nbonds[k];
            for (int l1 = 1; l1 <= nj; ++l1) {
                int l = ibonds[l1][j];
                if (l == i || mb[l] < 3) continue;
                for (int m1 = 1; m1 <= nk; ++m1) {
                    int m = ibonds[m1][k];
                    if (m == i || mb[m] < 3) continue;
                    int nm = nbonds[m];
                    bool found = false;
                    for (int n1 = 1; n1 <= nm; ++n1)
                        if (ibonds[n1][m] == l) { found = true; break; }
                    if (found) {
                        ir5[i] = ir5[j] = ir5[k] = ir5[l] = ir5[m] = ++nrings;
                        return;
                    }
                }
            }
        }
    }
}

// F90 749-803 -arom
// True if the II-JJ pi bond completes an aromatic ring, i.e. two other pi
// bonds already exist in that ring.
bool arom(int ii, int jj, const std::vector<int>& mpii) {
    int ni = nbonds[ii], nj = nbonds[jj];
    for (int i = 1; i <= ni; ++i) {
        int iia = ibonds[i][ii];
        if (iia == jj || mpii[iia] == 0) continue;
        int iib = mpii[iia];
        for (int j = 1; j <= nj; ++j) {
            int jja = ibonds[j][jj];
            if (jja == ii || mpii[jja] == 0) continue;
            int jjb = mpii[jja];
            int njb = nbonds[jjb];
            for (int j2 = 1; j2 <= njb; ++j2)
                if (ibonds[j2][jjb] == iib) return true;
        }
    }
    return false;
}

// F90 804-899 -arom2
// True if the II-JJ pi bond is the second pi bond of an aromatic ring; the
// remaining (third) bond is then assigned by rule 1 (arom).
bool arom2(int ii, int jj, const std::vector<int>& mpii) {
    int ni = nbonds[ii], nj = nbonds[jj];
    for (int i = 1; i <= ni; ++i) {
        int iia = ibonds[i][ii];
        if (iia == jj || mpii[iia] == 0) continue;
        int iib = mpii[iia];
        int nib = nbonds[iib];
        for (int j = 1; j <= nj; ++j) {
            int jja = ibonds[j][jj];
            if (jja == ii) continue;
            for (int i1 = 1; i1 <= nib; ++i1) {
                int iic = ibonds[i1][iib];
                for (int j1 = 1; j1 <= nbonds[jja]; ++j1)
                    if (iic == ibonds[j1][jja]) return true;
            }
        }
    }
    for (int j = 1; j <= nj; ++j) {
        int jja = ibonds[j][jj];
        if (jja == ii || mpii[jja] == 0) continue;
        int jjb = mpii[jja];
        int njb = nbonds[jjb];
        for (int i = 1; i <= ni; ++i) {
            int iia = ibonds[i][ii];
            if (iia == jj) continue;
            for (int j1 = 1; j1 <= njb; ++j1) {
                int jjc = ibonds[j1][jjb];
                for (int i1 = 1; i1 <= nbonds[iia]; ++i1)
                    if (jjc == ibonds[i1][iia]) return true;
            }
        }
    }
    return false;
}

// F90 1-618 -build the Lewis structure (sigma bonds, lone pairs, pi bonds).
void chklew(std::vector<int>& mb, std::vector<int>& numbon, int& l,
            int large, bool debug) {
    // large : number of atoms to dump in DEBUG output (20, or 10^6 for LARGE)
    // debug : TRUE when the GEOCHK keyword is present

    // F90 39-41 -(re)allocate the Lewis element table; LEWIS_TOT counts
    // entries.  (The Fortran allocates (2,norbs); rows 1-2 are all that is used
    // by this module and its callers.)
    if (!Lewis_elem.empty()) Lewis_elem.clear();
    Lewis_elem.assign(3, std::vector<int>(norbs + 1, 0));
    Lewis_tot = 0;

    bool big = keywrd.find(" LARGE") != std::string::npos;  // F90 42
    bool graphi = false;                                    // F90 43
    bool first_pi = true;                                   // F90 44
    int npi = 0;                                            // F90 30 / 230
    int not_aromatic_pi = 0;                                // F90 29 / 231

    // F90 45-49 -IPI(6*numat), IR5(numat), MPII(numat).  The Fortran checks
    // the allocation status and calls memory_error(); std::vector failure
    // would throw, which cannot be routed through that path.
    std::vector<int> ipi(6 * numat + 1, 0);
    std::vector<int> ir5(numat + 1, 0);
    std::vector<int> mpii(numat + 1, 0);

    // IPI is grown on demand instead of silently overrunning the Fortran
    // bound 6*numat (writing past a std::vector would be undefined behaviour).
    auto ipi_put = [&ipi](int idx, int value) {
        if (idx >= static_cast<int>(ipi.size())) ipi.resize(idx + 1, 0);
        ipi[idx] = value;
    };

    // F90 50-57 -orbitals and core electrons per atom.
    for (int i = 1; i <= numat; ++i) {
        ib[i] = nlast[i] - nfirst[i] + 1;
        iz[i] = static_cast<int>(std::lround(tore[nat[i]]));  // Nint(tore(nat(i)))
    }

    // F90 76-93 -explicit charges written in the atom labels.  Column 7 is
    // blanked while scanning columns 1-6 / 1-12, then restored.
    int j = 0;
    for (int i = 1; i <= natoms; ++i) {
        if (labels[i] != 99) ++j;
        std::string plus_seven = " ";
        const bool has_col7 = static_cast<int>(txtatm[i].size()) >= 7;
        if (has_col7) {
            plus_seven = txtatm[i].substr(6, 1);
            txtatm[i].replace(6, 1, " ");
        }
        if (txtatm[i].substr(0, 6).find('+') != std::string::npos) {
            ++icharges; iz[j]--; ions[j] = 1;    // F90 83-86
        } else if (txtatm[i].substr(0, 12).find('-') != std::string::npos) {
            ++icharges; iz[j]++; ions[j] = -1;   // F90 87-91
        }
        if (has_col7) txtatm[i].replace(6, 1, plus_seven);  // F90 92
    }

    // F90 97-110 -fill the "d" shell: NDELEC(I) electrons are not part of the
    // Lewis structure, so they are registered as lone pairs.
    for (int i = 1; i <= numat; ++i) {
        if (ib[i] != nlast[i] - nfirst[i] + 1) {
            int z = nat[i];
            if (ndelec[z] != 0)
                for (int k = 1; k <= ndelec[z] / 2; ++k)
                    add_Lewis_element(i, 0, 0, numbon[2]);
        }
    }

    // F90 112-123 -debug dump before the sigma framework.
    if (debug) {
        print_ld_text(" Arrays before making Sigma Framework (ATOM CORE ORBS BONDS)");
        int jd = std::min(numat, large);
        for (int k = 1; k <= jd; k += 20) {
            l = std::min(k + 19, jd);
            print_i4_row(k, l, [](int i) { return i; });
            print_i4_row(k, l, [&](int i) { return iz[i]; });
            print_i4_row(k, l, [&](int i) { return ib[i]; });
            print_i4_row(k, l, [&](int i) { return nbonds[i]; });
        }
    }
    numbon[1] = numbon[2] = numbon[3] = 0;  // F90 124 -numbon = 0

    // F90 128-144 -sigma framework.  A sigma bond is created only when the
    // second atom has the higher number (prevents double counting) and both
    // atoms still have an unused orbital.
    for (int ii = 1; ii <= numat; ++ii) {
        mb[ii] = 0;
        int nbii = nbonds[ii];
        for (int i = 1; i <= nbii; ++i) {
            int jj = ibonds[i][ii];
            if (jj >= ii && ib[ii] > 0 && ib[jj] > 0)
                add_Lewis_element(ii, jj, 0, numbon[1]);
        }
    }

    // F90 145-149 -label 1000: if every atom is saturated, L = 0 and return
    // (label 1100).  The tail below reproduces the L bookkeeping.
    {
        bool all_saturated = true;
        for (int i = 1; i <= numat; ++i)
            if (ib[i] != 0) { all_saturated = false; break; }
        if (all_saturated) goto finish;
    }

    // F90 151-161 -debug dump after the sigma framework.
    if (debug) {
        print_ld_text(" Core charges after making Sigma Framework (ATOM CORE ORBS)");
        int jd = std::min(numat, large);
        for (int k = 1; k <= jd; k += 20) {
            l = std::min(k + 19, jd);
            print_i4_row(k, l, [](int i) { return i; });
            print_i4_row(k, l, [&](int i) { return iz[i]; });
            print_i4_row(k, l, [&](int i) { return ib[i]; });
        }
    }

    // F90 162-183 -user-supplied pi bonds (keyword SETPI).
    // NOTE: SETPI input is not wired in this C++ port (readmo.cpp does not fill
    // PIBONDS, only PIBONDS_TXT); when PIBONDS is populated, row `loop` holds
    // (atom_i, atom_j) and a zero atom_i terminates the list, as in Fortran.
    if (keywrd.find(" SETPI") != std::string::npos) {
        for (int loop = 1; loop <= numat; ++loop) {
            if (static_cast<int>(pibonds.size()) <= loop ||
                static_cast<int>(pibonds[loop].size()) < 2)
                break;
            int ii = pibonds[loop][0];
            if (ii == 0) break;  // F90 168
            int jj = pibonds[loop][1];
            add_Lewis_element(ii, jj, 0, numbon[3]);
            if (ib[ii] < 0 || ib[jj] < 0) {
                if (keywrd.find(" RESEQ") == std::string::npos &&
                    keywrd.find(" SITE=") == std::string::npos &&
                    keywrd.find(" 0SCF") == std::string::npos) {
                    // F90 173-178 -dynamically sized I edit descriptors:
                    // width = 2 + floor(log10(atom_number*1.01)).
                    int w1 = 2 + static_cast<int>(std::log10(ii * 1.01));
                    int w2 = 2 + static_cast<int>(std::log10(jj * 1.01));
                    char buf[512];
                    std::snprintf(buf, sizeof buf,
                                  " Impossible user-supplied pi bond between %s%*d"
                                  " and %s%*d",
                                  elemnt[nat[ii]].c_str(), w1, ii,
                                  elemnt[nat[jj]].c_str(), w2, jj);
                    mopend(std::string(buf));
                    return;
                }
            }
        }
    }

    // F90 188-209 -assign all lone pairs.  Two cases:
    //   IZ > IB: normal lone pair (two electrons from the atom, one orbital);
    //   IZ < IB: "virtual" lone pair (empty orbital).
    l = 1;  // F90 188
    for (int i = 1; i <= numat; ++i) {
        if (iz[i] > ib[i] && ib[i] > 0) {
            int jl = std::min(ib[i], iz[i] - ib[i]);
            for (int k = 1; k <= jl; ++k)
                add_Lewis_element(i, 0, 0, numbon[2]);
        }
        if (iz[i] < ib[i] && ib[i] > 0) {
            int jl = ib[i] - iz[i];
            for (int k = 1; k <= jl; ++k)
                add_Lewis_element(0, i, 0, numbon[2]);  // element_type unused
        }
    }

    // F90 211-216 -label 1010: same early-out as 1000.
    {
        bool all_saturated = true;
        for (int i = 1; i <= numat; ++i)
            if (ib[i] != 0) { all_saturated = false; break; }
        if (all_saturated) goto finish;
    }

    // F90 217-226 -debug dump after the lone pairs.
    if (debug) {
        print_ld_text(" Core charges after adding Lone Pairs (ATOM CORE ORBS)");
        int jd = std::min(numat, large);
        for (int k = 1; k <= jd; k += 20) {
            l = std::min(k + 19, jd);
            print_i4_row(k, l, [](int i) { return i; });
            print_i4_row(k, l, [&](int i) { return iz[i]; });
            print_i4_row(k, l, [&](int i) { return ib[i]; });
        }
    }

    // F90 230-251 -third set: pi bonds attached to five-membered rings.
    npi = 0;
    not_aromatic_pi = 0;
    for (int ii = 1; ii <= numat; ++ii) {
        mpii[ii] = 0;
        ir5[ii] = 0;
        l = 0;
        if (ib[ii] != 0) {
            int nbii = nbonds[ii];
            for (int i = 1; i <= nbii; ++i)
                if (ib[ibonds[i][ii]] != 0) ++l;
        }
        mb[ii] = l;
    }
    for (int i = 1; i <= numat; ++i)
        if (mb[i] > 2 && ir5[i] == 0) ring5(i, mb, ir5);

    // F90 252-288 -open the five-membered rings that are fused to another
    // ring: create pi bonds between ring atoms and atoms outside the ring.
    for (int ii = 1; ii <= numat; ++ii) {
        if (ir5[ii] != 0) {
            int nbii = nbonds[ii];
            for (int i = 1; i <= nbii; ++i) {
                int jj = ibonds[i][ii];
                if (ir5[jj] != ir5[ii]) {
                    if (mb[jj] != 0 && ib[jj] != 0 && ib[ii] > 0) {
                        // F90 266 -ISET is set here but overwritten at F90 293.
                        --mb[ii]; --mb[jj];
                        mpii[ii] = jj; mpii[jj] = ii;
                        add_Lewis_element(ii, jj, 0, numbon[3]);
                        if (debug && big)
                            std::printf(" PI BOND BETWEEN ATOMS%5d AND%5d\n",
                                        ii, jj);
                        ++npi; ipi_put(npi, ii);
                        ++npi; ipi_put(npi, jj);
                    }
                }
            }
        }
    }

    // ======================================================================
    // F90 289-618 -outer DO loop.  Each pass:
    //   293-369  count pi neighbours (MB); unless GRAPHI, walk open-ended
    //            (degree-1) pi systems and add ene/yne pi bonds;
    //   375-379  if nothing was added: stop when every atom is saturated
    //            (goto 1060 -> finish with L = 0), else go to 1020 (debug);
    //   403-452  rebuild MB and the candidate atom list IPI;
    //   456-553  rule loop (Rule 1 = arom, 2 = arom2, 3 = any bond) with the
    //            `go to 1030` restart when a bond is added under loop /= 1;
    //   554-618  if rules added nothing: open one ring bond (label 1040), which
    //            restarts this whole pass.
    // The Fortran `exit` (F90 605) leaves the loop and falls into the tail.
    // ======================================================================
    for (;;) {
        bool iset = false;  // F90 293

        // F90 294-311 -number of pi neighbours of each atom (MB).
        for (int ii = 1; ii <= numat; ++ii) {
            l = 0;
            if (ib[ii] != 0) {
                int nbii = nbonds[ii];
                for (int i = 1; i <= nbii; ++i) {
                    int jj = ibonds[i][ii];
                    if (ib[jj] != 0) {
                        // Special case C-S-X: treat S as di-valent.
                        if ((nat[ii] == 6 && nat[jj] == 16 && nbonds[jj] > 1) ||
                            (nat[ii] == 16 && nbonds[ii] > 1 && nat[jj] == 6))
                            continue;
                        ++l;
                    }
                }
            }
            mb[ii] = l;
        }

        // F90 312-320 -debug dump of MB.
        if (debug) {
            print_ld_text(" Number of atoms pi bonding to each atom (MB)");
            int jd = std::min(numat, large);
            for (int jj2 = 1; jj2 <= jd; jj2 += 20) {
                int k = std::min(jj2 + 19, numat);
                print_i4_row(jj2, k, [](int i) { return i; });
                print_i4_row(jj2, k, [&](int i) { return mb[i]; });
            }
        }

        // F90 324-369 -open-ended (ene / yne) pi systems.  Loop 1 starts from
        // carbon ends, loop 2 from any end (distinguishes C-C-C-N from N-C-C-C).
        if (!graphi) {
            for (int loop = 1; loop <= 2; ++loop) {
                for (int i1 = 1; i1 <= numat; ++i1) {
                    if (mb[i1] != 1) continue;
                    if (!(loop != 1 || std::lround(tore[nat[i1]] - 4.0) == 0))
                        continue;
                    int ii = i1;
                    int nbii = nbonds[ii];
                    for (int i = 1; i <= nbii; ++i) {
                        int jj = ibonds[i][ii];
                        if (natorb[nat[ii]] == 9 || natorb[nat[jj]] == 9)
                            continue;
                        while (mb[jj] != 0 && ib[jj] != 0 && ib[ii] > 0) {
                            iset = true;
                            ++not_aromatic_pi;
                            --mb[ii]; --mb[jj];
                            mpii[ii] = jj; mpii[jj] = ii;
                            add_Lewis_element(ii, jj, 0, numbon[3]);
                            if (debug && big)
                                std::printf(" PI BOND BETWEEN ATOMS%5d AND%5d\n",
                                            ii, jj);
                            ++npi; ipi_put(npi, ii);
                            ++npi; ipi_put(npi, jj);
                            if (ib[ii] <= 0 || ib[jj] <= 0) break;  // F90 362
                        }
                    }
                }
            }
        }

        // F90 375-379 -nothing added?  Then either finish (label 1060) or
        // fall into the 1020 debug block and try to open a ring.
        if (iset) continue;  // F90: the whole if-block is skipped, DO repeats

        {
            bool all_saturated = true;
            for (int i = 1; i <= numat; ++i)
                if (ib[i] != 0) { all_saturated = false; break; }
            if (all_saturated) goto finish;  // F90 379 -go to 1060
        }

        // F90 380-398 -label 1020: debug dump.
        if (debug) {
            print_ld_text(" After removing open-ended pi systems");
            print_ld_text(" Number of bonds to each atom (MB, IB)");
            int jd = std::min(numat, large);
            for (int jj2 = 1; jj2 <= jd; jj2 += 20) {
                int k = std::min(jj2 + 19, numat);
                print_i4_row(jj2, k, [](int i) { return i; });
                print_i4_row(jj2, k, [&](int i) { return mb[i]; });
                print_i4_row(jj2, k, [&](int i) { return ib[i]; });
            }
            print_ld_text(" Core charges (ATOM CORE ORBS)");
            int jd2 = std::min(numat, large);
            for (int k = 1; k <= jd2; k += 20) {
                l = std::min(k + 19, jd2);
                print_i4_row(k, l, [](int i) { return i; });
                print_i4_row(k, l, [&](int i) { return iz[i]; });
                print_i4_row(k, l, [&](int i) { return ib[i]; });
            }
        }

        // F90 403-415 -rebuild MB: pi neighbours of the remaining pi atoms.
        for (int ii = 1; ii <= numat; ++ii) {
            l = 0;
            if (ib[ii] != 0) {
                int nbii = nbonds[ii];
                for (int i = 1; i <= nbii; ++i)
                    if (ib[ibonds[i][ii]] != 0) ++l;
            }
            mb[ii] = l;
        }

        // F90 423-426 -first call initialisation of NOT_AROMATIC_PI.
        if (first_pi) {
            not_aromatic_pi = 0;
            first_pi = false;
        }

        // F90 427-441 -extend IPI with the pi neighbours of the atoms already
        // in the list (LOOP1: skip a neighbour already listed).
        l = npi;
        for (int i = 1; i <= npi; ++i) {
            int ni_loc = nbonds[ipi[i]];
            for (int jj2 = 1; jj2 <= ni_loc; ++jj2) {
                int k = ibonds[jj2][ipi[i]];
                if (ib[k] == 0) continue;
                bool dup = false;
                for (int m = 1; m <= l; ++m)
                    if (ipi[m] == k) { dup = true; break; }
                if (dup) continue;
                ++l;
                ipi_put(l, k);
            }
        }
        npi = l;

        // F90 442-452 -keep only unique atoms that still have free valence.
        l = 0;
        for (int i = 1; i <= npi; ++i) {
            bool dup = false;
            for (int jj2 = 1; jj2 <= l; ++jj2)
                if (ipi[jj2] == ipi[i]) { dup = true; break; }
            if (dup) continue;
            if (ib[ipi[i]] != 0) {
                ++l;
                ipi_put(l, ipi[i]);
            }
        }
        npi = l;

        // F90 456-553 -pi-bond rules with restart.  `loop` starts at 1 and
        // advances to 2, 3; if a bond is added while loop /= 1 the original
        // does `go to 1030` (loop = 1) and rescans.  If loop reaches 4 the
        // rule phase is over and the ring phase follows.
        for (int loop = 1;;) {
            bool iset2 = false;  // F90 458
            bool restart = false;
            for (int i1 = npi; i1 >= 1 && !restart; --i1) {
                int ii = ipi[i1];
                if (mb[ii] == 0) continue;  // F90 461
                int nbii = nbonds[ii];
                for (int i = 1; i <= nbii; ++i) {
                    int jj = ibonds[i][ii];
                    if (natorb[nat[ii]] == 9 || natorb[nat[jj]] == 9) continue;
                    if (mb[jj] == 0) continue;                  // F90 466
                    if (!(ib[ii] > 0 && ib[jj] > 0)) continue;  // F90 467
                    // F90 483-487 -Rule 1 / Rule 2 / Rule 3.
                    bool make_bond;
                    if (loop == 1) make_bond = arom(ii, jj, mpii);
                    else if (loop == 2) make_bond = arom2(ii, jj, mpii);
                    else make_bond = true;
                    if (!make_bond) continue;
                    iset2 = true;  // F90 488
                    mpii[ii] = jj; mpii[jj] = ii;
                    add_Lewis_element(ii, jj, 0, numbon[3]);
                    // F90 495-497 -acetylenic bond needs a second pi bond.
                    if (ib[ii] == 1 && ib[jj] == 1)
                        add_Lewis_element(ii, jj, 0, numbon[3]);
                    // F90 503-523 -add JJ and the neighbours of II / JJ.
                    ++npi; ipi_put(npi, jj);
                    graphi = true;
                    l = npi;
                    int ni_loc = nbonds[ii];
                    for (int jj2 = 1; jj2 <= ni_loc; ++jj2) {
                        int k = ibonds[jj2][ii];
                        if (ib[k] != 0) { ++l; ipi_put(l, k); }
                    }
                    ni_loc = nbonds[jj];
                    for (int jj2 = 1; jj2 <= ni_loc; ++jj2) {
                        int k = ibonds[jj2][jj];
                        if (ib[k] != 0) { ++l; ipi_put(l, k); }
                    }
                    npi = l;
                    // F90 528-538 -LOOP3: drop duplicates / exhausted atoms.
                    l = 0;
                    for (int i2 = 1; i2 <= npi; ++i2) {
                        bool dup = false;
                        for (int jj2 = 1; jj2 <= l; ++jj2)
                            if (ipi[jj2] == ipi[i2]) { dup = true; break; }
                        if (dup) continue;
                        if (ib[ipi[i2]] != 0) { ++l; ipi_put(l, ipi[i2]); }
                    }
                    npi = l;
                }
                // F90 543 -a bond was made under loop /= 1: restart at 1030.
                if (iset2 && loop != 1) restart = true;
            }
            if (restart) {  // F90 1030
                loop = 1;
                continue;
            }
            ++loop;               // F90 546
            if (loop == 4) break; // F90 547-551 (otherwise CYCLE)
        }

        // F90 554-555 -ring phase.
        graphi = false;
        npi = 0;

        // F90 567-604 -open one ring pi bond.  The scan is ordered: loop 1
        // (atoms adjacent to a 3-pi-bond atom), loop 2 (2-pi-bond atoms), loop 3
        // (any remaining pi system).  Finding a bond is F90 `go to 1040`.
        bool ring_found = false;
        int ring_ii = 0, ring_jj = 0;
        for (int loop = 1; loop <= 3 && !ring_found; ++loop) {
            // F90 568-580 -rebuild MB.
            for (int ii = 1; ii <= numat; ++ii) {
                l = 0;
                if (ib[ii] != 0) {
                    int nbii = nbonds[ii];
                    for (int i = 1; i <= nbii; ++i)
                        if (ib[ibonds[i][ii]] != 0) ++l;
                }
                mb[ii] = l;
            }
            for (int ii = 1; ii <= numat && !ring_found; ++ii) {
                if (!(mb[ii] != 0 && (mb[ii] != 3 || loop >= 3))) continue;
                int nbii = nbonds[ii];
                bool lnext = false;
                for (int i = 1; i <= nbii; ++i)
                    if (mb[ibonds[i][ii]] == 3) lnext = true;
                lnext = lnext || loop > 2;
                for (int i = 1; i <= nbii; ++i) {
                    int jj = ibonds[i][ii];
                    if (!(lnext && mb[jj] != 0 && (mb[jj] != 3 || loop > 2)))
                        continue;
                    if (natorb[nat[ii]] == 9 || natorb[nat[jj]] == 9) continue;
                    if (ib[ii] > 0 && ib[jj] > 0) {  // F90 599 -go to 1040
                        ring_found = true;
                        ring_ii = ii;
                        ring_jj = jj;
                        break;
                    }
                }
            }
        }
        if (!ring_found) break;  // F90 605 -exit the outer DO loop

        // F90 1040-616 -add the ring pi bond, then repeat the outer DO loop.
        mpii[ring_ii] = ring_jj;
        mpii[ring_jj] = ring_ii;
        add_Lewis_element(ring_ii, ring_jj, 0, numbon[3]);
        if (debug && big)
            std::printf(" PI BOND BETWEEN ATOMS%5d AND%5d\n", ring_ii, ring_jj);
        npi = 1;
        ipi_put(npi, ring_ii);
        ++npi;
        ipi_put(npi, ring_jj);
    }

finish:
    // ------------------------------------------------------------------
    // F90 620-657 -tail of CHKLEW.
    //
    // F90 622-626: if no atom has an unused orbital left, the Lewis structure
    //              is complete; L is returned as 0.
    // F90 627-641: otherwise optionally dump MB / core charges / IB (DEBUG).
    // F90 642-643: success; L = 1.
    // F90 644-652: label 1060 -the "number of aromatic carbon atoms" block.
    //              It is only reached through the `goto 1060` at F90:379, the
    //              early exit taken when the pi framework is already complete.
    //              The count (I) is never used afterwards and L is explicitly
    //              reset to 0, which is exactly what the `else` branch below
    //              produces, so no separate code path is required.
    // F90 654:     deallocate(ipi, ir5, mpii) -these are std::vector locals
    //              here and need no explicit deallocation.
    // ------------------------------------------------------------------
    bool any_unused_orbital = false;
    for (int i = 1; i <= numat; ++i) {
        if (ib[i] > 0) {
            any_unused_orbital = true;
            break;
        }
    }
    if (any_unused_orbital) {
        if (debug) {
            // F90 627-641 -dump of the post-pi-system bookkeeping.
            print_ld_text(" After removing closed pi systems");
            print_ld_text(" Number of bonds to each atom");
            int k = std::min(20, numat);
            print_i4_row(1, k, [](int i) { return i; });
            print_i4_row(1, k, [&](int i) { return mb[i]; });
            print_ld_text(" Core charges (ATOM CORE ORBS)");
            int j = std::min(numat, large);
            for (int k1 = 1; k1 <= j; k1 += 20) {
                l = std::min(k1 + 19, j);  // F90 reuses L as the row bound
                print_i4_row(k1, l, [](int i) { return i; });
                print_i4_row(k1, l, [&](int i) { return iz[i]; });
                print_i4_row(k1, l, [&](int i) { return ib[i]; });
            }
        }
        l = 1;  // F90 642
    } else {
        l = 0;  // F90 625
    }

    // The variables below drive the control flow above but their final values
    // have no effect after CHKLEW returns.
    (void)big; (void)graphi; (void)first_pi; (void)not_aromatic_pi; (void)npi;
}
