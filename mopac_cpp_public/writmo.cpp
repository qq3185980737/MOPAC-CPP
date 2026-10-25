// writmo.cpp — C++ translation of "writmo.F90" (1469 lines; common output paths).
// WRITMO prints out most of the job results. It must not alter parameters.
//
// 【中文导读】writmo() —— 结果输出（对应 writmo.F90）
// ----------------------------------------------------------------------------
// 作用：把一次计算的主要结果打印到主输出（.out），函数不得修改任何参数。
// 覆盖范围（2026-10 复核）：主结果摘要（方法/版本/日期/SCF 与优化状态）、
// 压力/应力修正、FINAL HEAT OF FORMATION、DISP 能量分解、COSMO 默认腔体、
// 电离能/HOMO-LUMO、优化后几何与内坐标、Hessian、点群、本征值、原子电荷、
// 偶极矩、FOCK/DENS、AO 电子布居、PI 键级、UHF 自旋分析、键级、1ELE/ENPART、
// MULLIK/GRAPH 等常用表格均已移植。
// 尚缺（详见文件末尾 TODO 清单）：MOZYME 局域轨道输出、超精细耦合、pKa、
// PDB/图形输出、周期体系元胞/压力输出、公式与参数使用表等专用块。
//
// 单位约定：内部能量多用 eV；对外生成热/能量用 kcal/mol（fpc_9=23.060529
// kcal/mol per eV）；TOTAL ENERGY 一行按 EV 打印。输出宽度/格式刻意保持与
// 原版一致（%10s、%17.5f 等），以便与回归输出做文本对比。
//
// 阅读提示：函数开头用 flepo[]/iter_msg[] 查表打印优化/SCF 状态；之后按
// 关键字条件决定打印哪些段落，“输出内容随输入关键字变化”是正常现象。
// 原版 writmo.F90 共 1469 行，可作为对照逐段核对。
// ============================================================================

#include "writmo.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "cosmo.h"
#include "big_swap.h"          // l_control
#include "bonds.h"
#include "denrot.h"
#include "dimens.h"
#include "gmetry.h"
#include "write_cell.h"
#include "write_pressure.h"
#include "write_unit_cell_hof.h"
#include "volume.h"
#include "symtry.h"
#include "prtgra.h"
#include "geoutg.h"
#include "modchg.h"
#include "modgra.h"
#include "write_paras_used.h"
#include "dipole.h"
#include "dipole_for_MOZYME.h"
#include "dot.h"
#include "elemts_C.h"
#include "chrge.h"
#include "ijbo.h"
#include "enpart.h"
#include "meci_C.h"
#include "mecip.h"
#include "mullik.h"
#include "parameters_C.h"
#include "post_scf_corrections.h"
#include "funcon_C.h"
#include "maps_C.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "prtlmo.h"
#include "local_for_MOZYME.h"
#include "eigen.h"
#include "convert_storage.h"
#include "mpcsyb.h"
#include "timout.h"
#include "to_screen.h"
#include "vecprt.h"
#include "volume.h"
#include "wrttxt.h"
#include "symmetry_C.h"
#include "symtrz.h"
#include "diagg.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace cosmo_C;
using namespace maps_C;
using namespace molkst_C;

// Not-yet-ported helpers called by writmo.
void geout(int mode);
void PM7_TS();
void fdate(std::string& s);

// flepo(1..17): geometry-optimization outcome messages.
static const char* flepo[18] = {
    "",
    " 1SCF WAS SPECIFIED, SO BFGS WAS NOT USED",
    " GRADIENTS WERE INITIALLY ACCEPTABLY SMALL",
    " HERBERTS TEST WAS SATISFIED IN BFGS",
    " THE LINE MINIMIZATION FAILED TWICE IN A ROW.   TAKE CARE!",
    " BFGS FAILED DUE TO COUNTS EXCEEDED. TAKE CARE!",
    " PETERS TEST WAS SATISFIED IN BFGS OPTIMIZATION",
    " THIS MESSAGE SHOULD NEVER APPEAR, THERE IS A BUG IN MOPAC",
    " GRADIENT TEST NOT PASSED, BUT FURTHER WORK NOT JUSTIFIED",
    " A FAILURE HAS OCCURRED, TREAT RESULTS WITH CAUTION!!",
    " GEOMETRY OPTIMIZED USING NLLSQ. GRADIENT NORM MINIMIZED",
    " GEOMETRY OPTIMIZED USING POWSQ. GRADIENT NORM MINIMIZED",
    " CYCLES EXCEEDED, GRADIENT NOT FULLY MINIMIZED IN NLLSQ",
    " 1SCF RUN AFTER RESTART.  GEOMETRY MIGHT NOT BE OPTIMIZED",
    " HEAT OF FORMATION MINIMIZED IN ONE LINE SEARCH",
    " GEOMETRY OPTIMISED USING EIGENVECTOR FOLLOWING (EF).",
    " NO PARAMETERS MARKED FOR OPTIMIZATION, SO 1SCF WAS USED",
    " GEOMETRY OPTIMISED USING EIGENVECTOR FOLLOWING (TS)."};
static const char* iter_msg[3] = {"", " SCF FIELD WAS ACHIEVED",
                                  " ++++----**** FAILED TO ACHIEVE SCF. ****----++++"};

void writmo() {
    bool lprtgra = keywrd.find(" GRAD") != std::string::npos && nvar > 0;
    (void)lprtgra;
    int kchrge = 0;   // nint(total charge), used by the archive summary

    // 根据方法标志拼出方法名（AM1/PM3/...）；先保留前导空格再左裁剪。

    std::string caltyp = "   ";
    if (method_am1) caltyp = "       AM1";
    else if (method_pm3) caltyp = "       PM3";
    else if (method_mndo) caltyp = "      MNDO";
    else if (method_mndod) caltyp = "     MNDO-D";
    else if (method_rm1) caltyp = "       RM1";
    else if (method_pm6) caltyp = "       PM6";
    else if (method_pm7) caltyp = "       PM7";
    // left-trim caltyp
    size_t sp = caltyp.find_first_not_of(' ');
    if (sp != std::string::npos) caltyp = caltyp.substr(sp);

    std::printf("\n ----");
    for (int i = 0; i < 15; ++i) std::printf("-----");
    std::printf("\n");
    wrttxt(iw);

    double gnorm = 0.0;
    if (nvar != 0) gnorm = std::sqrt(dot(grad, grad, nvar));
    (void)gnorm;

    if (iflepo == 15 && keywrd.find(" TS") != std::string::npos) iflepo = 17;
    if (iflepo == 0) iflepo = 7;
    if (iflepo == 1 && nvar == 0 && keywrd.find("1SCF") == std::string::npos)
        iflepo = 16;
    std::printf("\n\n    %s\n", flepo[iflepo]);
    if (iscf < 1) iscf = 1;
    std::printf("    %s\n", iter_msg[iscf]);
    std::printf("\n\n%30s%s CALCULATION\n", "", caltyp.c_str());
    std::printf("%55sMOPAC2016 (Version: %s)\n", "", verson.c_str());

    std::string idate;
    fdate(idate);
    std::printf("%55s%s\n", "", idate.c_str());

    if (iscf == 2) {
        std::printf(" \n \n FOR SOME REASON THE SCF CALCULATION FAILED.\n \n");
        std::printf(" THE RESULTS WOULD BE MEANINGLESS, SO WILL NOT BE PRINTED.\n");
        std::printf(" TRY TO FIND THE REASON FOR THE FAILURE BY USING 'PL'.\n \n");
        std::printf(" CHECK YOUR GEOMETRY AND ALSO TRY USING SHIFT OR PULAY. \n");
        geout(1);
        return;
    }

    // External pressure contribution.
    if (std::abs(pressure) > 1e-4) {
        if (id == 1) {
            double ln = 0.0;
            for (int d = 0; d < 3; ++d) ln += tvec[d][0] * tvec[d][0];
            escf = escf + pressure * std::sqrt(ln);
        } else if (id == 3) {
            double tv[10];
            for (int d = 0; d < 3; ++d)
                for (int j = 0; j < 3; ++j) tv[d * 3 + j] = tvec[d][j];
            escf = escf + pressure * volume(tv, 3);
        }
    }

    if (use_ref_geo) {
        std::printf("\n\n%10sFINAL H.O.F PLUS STRESS =%17.5f KCAL/MOL  =%14.5f KJ/MOL\n",
                    "", escf, escf * 4.184);
        stress = stress * density;
        std::printf("%10sFINAL STRESS            =%17.5f KCAL/MOL  =%14.5f KJ/MOL\n",
                    "", stress, stress * 4.184);
        std::printf("%10sFINAL HEAT OF FORMATION =%17.5f KCAL/MOL  =%14.5f KJ/MOL\n",
                    "", escf - stress, (escf - stress) * 4.184);
        double sum = 0.0, rms = 0.0;
        for (int i = 1; i <= numat; ++i) {
            double dx = geo[0][i] - geoa[0][i];
            double dy = geo[1][i] - geoa[1][i];
            double dz = geo[2][i] - geoa[2][i];
            sum += std::sqrt(dx * dx + dy * dy + dz * dz);
            rms += dx * dx + dy * dy + dz * dz;
        }
        double distortion = sum / numat;
        std::printf("%10sTOTAL DISTORTION        =%17.5f Angstroms\n", "", sum);
        std::printf("%10sAVERAGE DISTORTION      =%17.5f Angstroms per atom (all atoms)\n",
                    "", distortion);
        std::printf("%10sRMS DISTORTION          =%17.5f Angstroms per atom (all atoms)\n",
                    "", std::sqrt(rms / numat));
    } else {
        if (keywrd.find(" PM7-TS") != std::string::npos) {
            PM7_TS();
            return;
        }
        stress = -1.0;
        std::printf("\n\n\n\n%10sFINAL HEAT OF FORMATION =%17.5f KCAL/MOL =%14.5f KJ/MOL\n",
                    "", escf, escf * 4.184);
    }

    // DISP：打印色散/氢键等能量分解，并给出加上校正项后的总生成热。

    if (keywrd.find(" DISP") != std::string::npos) {
        // Full DISP energy decomposition (F90 writmo L205-230).
        std::printf("\n%10sTOTAL ENERGY            =%17.5f KCAL/MOL\n", "",
                    (elect + enuclr + solv_energy) * funcon_C::fpc_9);
        std::printf("%10sENERGY OF ATOMS         =%17.5f KCAL/MOL\n", "", atheat);
        std::printf("%10s                    SUM =%17.5f KCAL/MOL\n", "",
                    (elect + enuclr) * funcon_C::fpc_9 + atheat +
                        solv_energy * funcon_C::fpc_9);
        if (std::abs(hpress) > 1.0e-5)
            std::printf("%10sENERGY DUE TO PRESSURE  =%17.5f KCAL/MOL\n", "", hpress);
        std::printf("%10sDISPERSION ENERGY       =%17.5f KCAL/MOL\n", "", E_disp);
        if (E_hb < -1.0e-5)
            std::printf("%10sH-BOND ENERGY           =%17.5f KCAL/MOL\n", "", E_hb);
        if (std::abs(nsp2_corr) > 1.0e-5)
            std::printf("%10sMM CORRECTION FOR >N-   =%17.5f KCAL/MOL\n", "", nsp2_corr);
        if (std::abs(Si_O_H_corr) > 1.0e-5)
            std::printf("%10sMM CORR. FOR Si-O-H     =%17.5f KCAL/MOL\n", "", Si_O_H_corr);
        if (std::abs(sum_dihed) > 1.0e-5)
            std::printf("%10sMM CORR. FOR -CO-NH-    =%17.5f KCAL/MOL\n", "", sum_dihed);
        double dsum = (elect + enuclr) * funcon_C::fpc_9 + atheat + hpress +
                      solv_energy * funcon_C::fpc_9 + nsp2_corr +
                      Si_O_H_corr + sum_dihed + E_disp + E_hb;
        std::printf("%30sSUM =%17.5f KCAL/MOL\n\n", "", dsum);
        if (std::abs(dsum - escf) > 1.0e-2)
            std::printf(" WARNING - An energy term is missing!\n");
        if (N_Hbonds > 0)
            std::printf("%10sNo. OF HYDROGEN BONDS   =%11d%7s\n", "", N_Hbonds,
                        "(H-bond Energy < -1.0 Kcal/mol)");
        if (E_hh > 1.0e-5)
            std::printf("%10sH - H CORRECTION ENERGY =%17.5f KCAL/MOL\n", "", E_hh);
        if (keywrd.find(" DISP(") != std::string::npos) {
            l_control("PRT", 3, 1);
            double corr = 0.0;
            post_scf_corrections(corr, false);
        }
    }
    // F90 L231-239: large-HOF warning and screen messages.
    if (numat > 1 && iscf == 1 && escf > 1.e4 &&
        keywrd.find(" CHECK") == std::string::npos)
        std::printf("\n\n%10sCalculated Heat of Formation is very large, re-run using keyword 'CHECK'\n\n", "");
    to_screen(" Job: " + jobnam);
    {
        char buf[160];
        std::snprintf(buf, sizeof(buf),
                      "Final heat of formation = %16.5f kcal/mol", escf);
        to_screen(buf);
        std::snprintf(buf, sizeof(buf),
                      "GRADIENT NORM           =%17.5f", gnorm);
        to_screen(buf);
    }

    // F90 writmo.F90:233: heat of formation per unit cell (periodic).
    if (id == 3) write_unit_cell_HOF(iw);
    if (keywrd.find(" EPS") != std::string::npos)
        std::printf("%10sVAN DER WAALS AREA      =%14.2f SQUARE ANGSTROMS\n", "", area);

    // F90 writmo.F90:243: STATE is only appended for open-shell / non-integral
    // occupations; closed-shell RHF prints the plain TOTAL ENERGY line.
    const bool has_state = (uhf || nopen != nclose || fract < 1.99) &&
                           symmetry_C::state_Irred_Rep.find_first_not_of(' ') != std::string::npos;
    if (has_state)
        std::printf("\n\n%10sTOTAL ENERGY            =%17.5f EV   STATE:  %2d %s %s\n",
                    "", elect + enuclr + solv_energy, symmetry_C::state_QN,
                    symmetry_C::state_spin.c_str(), symmetry_C::state_Irred_Rep.c_str());
    else
        std::printf("\n\n%10sTOTAL ENERGY            =%17.5f EV\n", "",
                    elect + enuclr + solv_energy);
    // F90 writmo.F90:254：先确定点群并对称化轨道/本征值，再打印点群。
    if (norbs > 0 && !mozyme && numat < 200) {
        std::vector<double> cflat((size_t)norbs * norbs + 1, 0.0);
        for (int col = 1; col <= norbs; ++col)
            for (int row = 1; row <= norbs; ++row)
                cflat[(col - 1) * norbs + row] = c[row][col];
        symtrz(&cflat[1], &eigs[1], 1, 1);
        if (moperr) return;
        for (int col = 1; col <= norbs; ++col)
            for (int row = 1; row <= norbs; ++row)
                c[row][col] = cflat[(col - 1) * norbs + row];
    }
    if (id == 0 && !mozyme)
        std::printf("%10sELECTRONIC ENERGY       =%17.5f EV  POINT GROUP:     %-4s\n",
                    "", elect, symmetry_C::name.c_str());
    else
        std::printf("%10sELECTRONIC ENERGY       =%17.5f EV\n", "", elect);
    std::printf("%10sCORE-CORE REPULSION     =%17.5f EV\n", "", enuclr);
    if (iseps)
        std::printf("%10sDIELECTRIC ENERGY       =%17.5f EV\n", "", ediel);

    // COSMO 默认腔体：即使用户没写 EPS，也对有限体系（非 MOZYME、原子数<700）
    // 用默认介电常数计算并打印面积/体积。

    // F90 writmo.F90: default COSMO cavity for finite, non-polymeric systems:
    // even without an explicit EPS keyword the area/volume are computed and
    // printed (fepsi = (78.4-1)/(78.4+0.5), cosini + coscav).
    if (keywrd.find(" EPS") == std::string::npos && id == 0 &&
        !mozyme && numat < 700) {
        const double sum = 78.4;
        fepsi = (sum - 1.0) / (sum + 0.5);
        cosini(false);
        if (moperr) {
            moperr = false;
            area = 0.0;
            cosvol = 0.0;
        } else {
            coscav();
        }
    }
    if (moperr) {
        moperr = false;
    } else if (fepsi > 1e-10 && area > 1e-3) {
        std::printf("%10sCOSMO AREA              =%14.2f SQUARE ANGSTROMS\n", "", area);
        std::printf("%10sCOSMO VOLUME            =%14.2f CUBIC ANGSTROMS\n", "", cosvol);
    }
    if (lprtgra || gnorm > 1.e-15)
        std::printf("%10sGRADIENT NORM           =%17.5f\n", "", gnorm);
    // F90 L345-369: stationary-point check + NOANCI hint.
    bool still = true;
    if (latom == 0) {
        if (keywrd.find(" AIDER") == std::string::npos) {
            if (keywrd.find(" 1SCF") == std::string::npos ||
                keywrd.find(" GRAD") != std::string::npos) {
                double gsum = std::sqrt(dot(dxyz, dxyz, 3 * numat));
                if (gsum > std::max(5.0, 5 * gnorm * gnorm) &&
                    gnorm < 2.0 && nclose == nopen && id == 0 &&
                    keywrd.find("NOANCI") == std::string::npos) {
                    if (nvar != 1 ||
                        keywrd.find(" GRAD") != std::string::npos ||
                        keywrd.find("DERIV") != std::string::npos) {
                        std::printf("%9s WARNING -- GEOMETRY IS NOT AT A STATIONARY POINT\n", "");
                        still = false;
                    }
                }
            }
        }
    }
    // Ionization potential and HOMO/LUMO energies.

    // 中文：电离能用 Koopmans 近似（-HOMO）；开壳层取 alpha/beta 中较高者，
    // 并对开壳层 doublet 做 rjkab 修正（见下方 eionis 计算）。

    double eionis = 0.0;
    if (nalpha > 0 && nbeta > 0)
        eionis = -std::max(eigs[nalpha], eigb[nbeta]);
    else if (nelecs == 1)
        eionis = -eigs[1];
    else if (nelecs > 1) {
        if (nclose > 0) eionis = -eigs[nclose];
        if (nopen > 0) eionis = std::min(eionis, -eigs[nopen]);
    }
    {
        int iop = nclose;
        if (fract > 1.99) iop = nopen;
        int nopn = nopen - iop;
        if (!mozyme) {
            if (nopn == 1 && meci_C::rjkab.size() > 1)
                eionis += 0.5 * meci_C::rjkab[1][1];  // doublet IP correction
            if (std::abs(eionis) > 1e-5 && nopn < 2)
                std::printf("\n%10sIONIZATION POTENTIAL    =  %16.6f EV\n", "", eionis);
            if (uhf) {
                if (nalpha >= 1)
                    std::printf("%10sALPHA SOMO LUMO (EV)    =  %13.3f%7.3f\n", "",
                                eigs[nalpha],
                                (nalpha + 1 <= norbs) ? eigs[nalpha + 1] : 0.0);
                if (nbeta >= 1)
                    std::printf("%10sBETA  SOMO LUMO (EV)    =  %13.3f%7.3f\n", "",
                                eigb[nbeta],
                                (nbeta + 1 <= norbs) ? eigb[nbeta + 1] : 0.0);
            } else if (nopen == nclose) {
                if (nopen >= 1)
                    std::printf("%10sHOMO LUMO ENERGIES (EV) =  %13.3f%7.3f\n", "",
                                eigs[nopen],
                                (nopen + 1 <= norbs) ? eigs[nopen + 1] : 0.0);
            } else if (nopn == 1) {
                if (nclose >= 1)
                    std::printf("%10sHOMO (SOMO) LUMO (EV)   =%13.3f (%7.3f)%7.3f\n",
                                "", eigs[nclose],
                                (nclose + 1 <= norbs) ? eigs[nclose + 1] : 0.0,
                                (nclose + 2 <= norbs) ? eigs[nclose + 2] : 0.0);
                else
                    std::printf("%10s     (SOMO) LUMO (EV)   =      (%6.3f)%7.3f\n",
                                "",
                                (nclose + 1 <= norbs) ? eigs[nclose + 1] : 0.0,
                                (nclose + 2 <= norbs) ? eigs[nclose + 2] : 0.0);
            }
        }
        if (uhf) {
            std::printf("%10sNO. OF ALPHA ELECTRONS  =%11d\n", "",
                        nalpha + (int)((nalpha_open - nalpha) * fract + 0.5));
            std::printf("%10sNO. OF BETA  ELECTRONS  =%11d\n", "",
                        nbeta + (int)((nbeta_open - nbeta) * fract + 0.5));
        } else {
            if (mozyme) std::printf("\n");
            std::printf("%10sNO. OF FILLED LEVELS    =%11d\n", "", nopen - nopn);
            if (nopn != 0)
                std::printf("%10sAND NO. OF OPEN LEVELS  =%11d\n", "", nopn);
        }
    }
    if (mol_weight > 0.1)
        std::printf("%10sMOLECULAR WEIGHT        =%16.4f\n", "", mol_weight);
    // F90 writmo.F90:447: re-derive Cartesian coord from the internal geometry
    // right before dimens/printing (symmetry handling may have rotated coord).
    gmetry(geo, coord);
    dimens(coord, iw);   // F90: call dimens(coord, iw) after MOLECULAR WEIGHT
    // F90 writmo.F90:449-458: periodic systems print unit-cell volume/density
    // and the cell/pressure tables.
    if (id == 3) {
        double vecs[9];
        for (int ii = 1; ii <= 3; ++ii)
            for (int jj = 1; jj <= 3; ++jj)
                vecs[(jj - 1) * 3 + (ii - 1)] = tvec[ii][jj];
        const double vol = volume(vecs, 3);
        const double density = mol_weight * 1.0e24 / funcon_C::fpc_10 / vol;
        std::printf("\n%10sVOLUME OF UNIT CELL%17.3f CUBIC ANGSTROMS\n",
                    "", vol);
        std::printf("%10sDENSITY        %21.3f GRAMS/CC\n", "", density);
        write_cell(iw);
        write_pressure(iw);
    }
    std::printf("\n\n%10sSCF CALCULATIONS        =   %8d\n", "", nscf);
    // F90 writmo.F90:465: write(iw,*) after the SCF count (single blank record).
    std::printf(" \n");
    timout(iw);   // F90: call timout(iw) after SCF CALCULATIONS
    if (ndep != 0) symtry();   // F90 writmo.F90:467
    // Sync optimized geometry back into xparam.
    for (int i = 1; i <= nvar; ++i)
        xparam[i] = geo[loc[2][i]][loc[1][i]];

    // FINAL POINT AND DERIVATIVES table.
    if (prt_gradients && (lprtgra || gnorm > 2.0)) {
        std::printf("\n\n\n%7sFINAL  POINT  AND  DERIVATIVES\n\n", "");
        if (mozyme) {
            prtgra();   // F90 writmo.F90:475 MOZYME gradient table
            if (maxtxt == 26) modgra();   // F90 writmo.F90:476 residue gradients
        } else {
            std::printf("   PARAMETER     ATOM    TYPE            VALUE       GRADIENT\n");
            const double degree = 57.29577951308232;
            for (int i = 1; i <= nvar; ++i) {
                int j = loc[2][i], k = loc[1][i];
                int l = labels[k];
                double xi = xparam[i];
                if (j != 1 && na[k] > 0) xi *= degree;
                const char* gtype = (j == 1 || na[k] == 0) ? "KCAL/ANGSTROM"
                                                          : "KCAL/RADIAN  ";
                const char* type;
                if (na[k] == 0)
                    type = (j == 1) ? "CARTESIAN X" : (j == 2) ? "CARTESIAN Y"
                                                                  : "CARTESIAN Z";
                else
                    type = (j == 1) ? "BOND       " : (j == 2) ? "ANGLE      "
                                                                : "DIHEDRAL   ";
                std::printf("%7d%11d  %s   %-11s%13.6f%13.6f  %s\n",
                            i, k, elemts_C::elemnt[l].c_str(), type, xi, grad[i], gtype);
            }
        }
    }

    // F90 writmo.F90:506-518: blank lines, geout(1) if prt_coords, then the
    // (i4,3x,A2,3x,3F16.9) Cartesian table.
    std::printf("\n\n\n\n");
    if (prt_coords) geout(1);
    if (prt_cart) {
        std::printf("\n%28s CARTESIAN COORDINATES\n\n", "");
        for (int i = 1; i <= numat; ++i) {
            std::string sym = elemts_C::elemnt[nat[i]];
            while (!sym.empty() && sym.back() == ' ') sym.pop_back();
            std::printf("%4d   %2s   %16.9f%16.9f%16.9f\n", i, sym.c_str(),
                        coord[0][i], coord[1][i], coord[2][i]);
        }
    }
    // PRTINT: interatomic distances (F90 L527-541).
    if (!mozyme && keywrd.find("PRTINT") != std::string::npos) {
        int npr = numat * (numat + 1) / 2;
        std::vector<double> rxyz(npr, 0.0);
        int l = -1;
        for (int i = 1; i <= numat; ++i)
            for (int j = 1; j <= i; ++j) {
                ++l;
                double dx = coord[0][i] - coord[0][j];
                double dy = coord[1][i] - coord[1][j];
                double dz = coord[2][i] - coord[2][j];
                rxyz[l] = std::sqrt(dx * dx + dy * dy + dz * dz);
            }
        std::printf("\n\n%10s  INTERATOMIC DISTANCES\n", "");
        vecprt(&rxyz[0], numat);
    }

    // Clip wild MO eigenvalues.
    if (!mozyme) {
        for (int i = 1; i <= norbs; ++i) {
            if (!(eigs[i] >= -999.0 && eigs[i] <= 1000.0)) eigs[i] = 1e-12;
            if (eigb[i] < -999.0 || eigb[i] > 1000.0) eigb[i] = 1e-12;
        }
    }
    // F90 writmo.F90:549: (/,a,/)
    if (!formula.empty()) std::printf("\n\n%s\n\n", formula.c_str());

    // Hessian matrix (upper triangle).
    if (keywrd.find(" HESSIAN") != std::string::npos) {
        std::printf("    Hessian matrix from geometry optimization\n");
        for (int i = 1; i <= nvar; ++i) {
            for (int j = 1; j <= i; ++j)
                std::printf("%10.4f", hesinv[(i - 1) * nvar + j - 1]);
            std::printf("\n");
        }
    }
    // MOLECULAR POINT GROUP.
    if (id == 0 && !mozyme)
        std::printf("\n\n      MOLECULAR POINT GROUP   :   %-4s\n",
                    symmetry_C::name.c_str());

    // MOZYME LMO output (F90 writmo.F90 L556-575).
    if (norbs * nelecs > 0 && mozyme) {
        const bool vec = keywrd.find(" VEC") != std::string::npos ||
                         keywrd.find(" ALLVEC") != std::string::npos;
        if (keywrd.find(" RE-LOC") != std::string::npos) {
            std::printf("  LMOs being Re-Localized\n\n");
            local_for_MOZYME("OCCUPIED");
            local_for_MOZYME("VIRTUAL");
        }
        if (vec && nelecs != 0) {
            if (keywrd.find(" EIGEN") != std::string::npos) {
                eigen(false, true);
            } else {
                // Official prtlmo prints eigs(1:noccupied) as left by the SCF:
                // those are the canonical/MOZYME eigenvalues of the final
                // diagonalisation, not a fresh Fock expectation value.
                prtlmo();
            }
        }
    }

    // EIGENVALUES list (unless eigenvectors requested via VEC/ALLVEC).
    // F90 writmo.F90 L819: the MO eigenvalue block is for non-MOZYME only;
    // MOZYME writes LMOs instead (see block above).
    if (norbs > 0 && !mozyme) {
        if (!(keywrd.find(" VEC") != std::string::npos ||
              keywrd.find(" ALLVEC") != std::string::npos)) {
            if (uhf) std::printf("\n\n%10sALPHA EIGENVALUES\n", "");
            else     std::printf("\n\n%18sEIGENVALUES\n", "");
            for (int i = 1; i <= norbs; ++i) {
                std::printf("%10.5f", eigs[i]);
                if (i % 8 == 0) std::printf("\n");
            }
            if (norbs % 8 != 0) std::printf("\n");
            if (uhf) {
                std::printf("\n\n%10s BETA EIGENVALUES\n", "");
                for (int i = 1; i <= norbs; ++i) {
                    std::printf("%10.5f", eigb[i]);
                    if (i % 8 == 0) std::printf("\n");
                }
                if (norbs % 8 != 0) std::printf("\n");
            }
        }
    }

    // NET ATOMIC CHARGES.
    // 中文：打印 Mulliken 净电荷；打印前会按需修正密度矩阵（见下方注释）。

    if (nelecs != 0) {
        // Correct density matrix, if necessary (F90 writmo.F90 L605-610).
        if ((nclose != nopen && std::fabs(fract - 2.0) > 1e-20 && fract > 1e-20) ||
            keywrd.find(" C.I.") != std::string::npos)
            mecip();
        if (prt_charges) {
            std::printf("\n\n%13s NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS\n\n", "");
            // F90 L614-625: table header with s-/p-/d-Pop labels.
            int nb = 0;
            for (int j = 1; j <= numat; ++j)
                nb = std::max(nb, parameters_C::natorb[nat[j]]);
            nb = (int)std::lround(std::sqrt((double)nb)) * 12;
            const char* pops = "s-Pop       p-Pop       d-Pop";
            if (maxtxt == 26)
                std::printf("  ATOM NO.                 TYPE                        CHARGE      No. of ELECS.   %s\n",
                            std::string(pops).substr(0, std::min(nb, 36)).c_str());
            else
                std::printf("  ATOM NO.   TYPE          CHARGE      No. of ELECS.   %s\n",
                            std::string(pops).substr(0, std::min(nb, 36)).c_str());
        }
        std::vector<double> q2(numat + 1, 0.0);
        chrge(p, q2);
        double sumq = 0.0;
        std::vector<double> rxyz(4, 0.0);
        for (int i = 1; i <= numat; ++i) {
            int l = nat[i];
            q[i] = parameters_C::tore[l] - q2[i];
            sumq += q[i];
            if (!l_atom[i]) continue;
            // F90 writmo.F90 L633-658: s-/p-/d-populations (rxyz(1:j)).
            int j = (int)std::lround(std::sqrt(parameters_C::natorb[l] + 1e-4));
            int k = nfirst[i];
            rxyz[1] = rxyz[2] = rxyz[3] = 0.0;
            if (j > 0) {
                if (mozyme) {
                    k = ijbo(i, i);
                    rxyz[1] = p[k + 1];
                    if (j > 1) {
                        rxyz[2] = p[k + 3] + p[k + 6] + p[k + 10];
                        if (j > 2)
                            rxyz[3] = p[k + 15] + p[k + 21] + p[k + 28] + p[k + 36] + p[k + 45];
                    }
                } else {
                    rxyz[1] = p[(k * (k + 1)) / 2];
                    if (j > 1) {
                        rxyz[2] = p[((k + 1) * (k + 2)) / 2] + p[((k + 2) * (k + 3)) / 2]
                                + p[((k + 3) * (k + 4)) / 2];
                        if (j > 2)
                            rxyz[3] = p[((k + 4) * (k + 5)) / 2] + p[((k + 5) * (k + 6)) / 2]
                                    + p[((k + 6) * (k + 7)) / 2] + p[((k + 7) * (k + 8)) / 2]
                                    + p[((k + 8) * (k + 9)) / 2];
                    }
                }
            }
            if (prt_charges) {
                // F90 writmo.F90:661-667:
                //   (I5,9X,A2,4X,F15.6,F14.4,3F12.5)   (or A30 label variant)
                std::string sym = elemts_C::elemnt[l];
                while (!sym.empty() && sym.back() == ' ') sym.pop_back();
                if (maxtxt == 26) {
                    std::string lab = elemts_C::elemnt[l] + "(" + txtatm[i] + ")";
                    std::printf("%5d         %-30s    %15.6f%14.4f", i, lab.c_str(),
                                q[i], q2[i]);
                } else {
                    std::printf("%5d         %2s    %15.6f%14.4f", i, sym.c_str(),
                                q[i], q2[i]);
                }
                for (int m = 1; m <= j; ++m) std::printf("%12.5f", rxyz[m]);
                std::printf("\n");
            }
        }
        if (maxtxt == 26) modchg();   // F90 writmo.F90:670 residue charge table
        kchrge = (int)std::lround(sumq);
        (void)kchrge;
    }
    // DIPOLE section (F90: dipole(p, coord, dumy, 1) after NET ATOMIC CHARGES).

    // 中文：偶极矩由 dipole.cpp 计算，依赖密度矩阵 p 与几何 coord。

    if (id == 0) {
        std::vector<double> dipvec(4, 0.0);
        if (mozyme) dipole_for_MOZYME(dipvec, 1);
        else dipole(p, coord, dipvec, 1);
    }
    // FOCK / DENS matrices (F90 L683-698).
    if (norbs > 0) {
        if (keywrd.find(" FOCK") != std::string::npos) {
            std::printf(" FOCK MATRIX \n");
            vecprt(&f[0], norbs);
        }
        if (nelecs != 0 && keywrd.find(" DENS") != std::string::npos) {
            std::printf("\n\n%20s DENSITY MATRIX IS \n", "");
            vecprt(&p[0], norbs);
        }
    }
    // ATOMIC ORBITAL ELECTRON POPULATIONS.
    if (nelecs != 0 && prt_pops) {
        std::printf("\n\n%10sATOMIC ORBITAL ELECTRON POPULATIONS\n\n", "");
        int maxorb = 0;
        for (int i = 1; i <= numat; ++i)
            maxorb = std::max(maxorb, nlast[i] - nfirst[i]);
        int ihead = (maxtxt == 0) ? 2 : (maxtxt == 26 ? 16 : maxtxt / 2 + 2);
        std::string hsp((size_t)ihead, ' ');
        if (maxorb == 8)
            std::printf("%s   Atom%s  s        px        py        pz      x^2-y^2     xz        z^2       yz        xy\n",
                        hsp.c_str(), hsp.c_str());
        else if (maxorb == 3)
            std::printf("%s   Atom%s  s        px        py        pz   \n",
                        hsp.c_str(), hsp.c_str());
        else
            std::printf("%s   Atom%s  s   \n", hsp.c_str(), hsp.c_str());
        if (mozyme) {
            // F90 writmo.F90 L724-743: pdiag(l) = p(ijbo(i,i)+k), then print by
            // global AO indices nfirst..nlast.
            std::vector<double> pdiag(norbs + 1, 0.0);
            int l = 0;
            for (int i = 1; i <= numat; ++i) {
                int j = ijbo(i, i);
                for (int k = 1; k <= MOZYME_C::iorbs[i]; ++k) {
                    ++l;
                    j += k;
                    pdiag[l] = p[j];
                }
            }
            for (int i = 1; i <= numat; ++i) {
                std::string sym = elemts_C::elemnt[nat[i]];
                while (!sym.empty() && sym.back() == ' ') sym.pop_back();
                std::printf("%5d %2s", i, sym.c_str());
                for (int j = nfirst[i]; j <= nlast[i]; ++j)
                    std::printf("%10.5f", pdiag[j]);
                std::printf("\n");
            }
        } else {
            for (int i = 1; i <= numat; ++i) {
                std::string sym = elemts_C::elemnt[nat[i]];
                while (!sym.empty() && sym.back() == ' ') sym.pop_back();
                std::printf("%5d %2s", i, sym.c_str());
                for (int j = nfirst[i]; j <= nlast[i]; ++j)
                    std::printf("%10.5f", p[(j * (j + 1)) / 2]);
                std::printf("\n");
            }
        }
    }
    // PI bond orders (F90 L758-761).
    if (keywrd.find(" PI") != std::string::npos) {
        std::printf("\n\n%10sSIGMA-PI BOND-ORDER MATRIX\n", "");
        denrot();
    }
    // UHF spin analysis: SZ and <S**2>.
    if (nelecs != 0 && uhf) {
        double na_u = nalpha + (nalpha_open - nalpha) * fract;
        double nb_u = nbeta + (nbeta_open - nbeta) * fract;
        sz = (na_u - nb_u) * 0.5;
        ss2 = sz * sz;
        int l = 0;
        for (int i = 1; i <= norbs; ++i) {
            for (int j = 1; j <= i; ++j) {
                ++l;
                pa[l] = pa[l] - pb[l];
                ss2 = ss2 + pa[l] * pa[l];
            }
            ss2 = ss2 - 0.5 * pa[l] * pa[l];
        }
        std::printf("\n\n%20s(SZ)    =%12.6f\n", "", sz);
        if (fract < 1e-5 || fract > 0.9999)
            std::printf("%20s(S**2)  =%12.6f\n", "", ss2);
        else
            std::printf("%10sAverage over configurations used, so (S**2) is not meaningful\n", "");
        // SPIN density matrix (F90 L781-784).
        if (keywrd.find(" SPIN") != std::string::npos) {
            std::printf("\n\n%10sSPIN DENSITY MATRIX\n", "");
            vecprt(&pa[0], norbs);
        }

        // AO spin populations.
        std::printf("\n\n%10sATOMIC ORBITAL SPIN POPULATIONS\n\n", "");
        int maxorb = 0;
        for (int i = 1; i <= numat; ++i)
            maxorb = std::max(maxorb, nlast[i] - nfirst[i]);
        if (maxorb == 8)
            std::printf("    Atom   Total    s     px     py     pz   x^2-y^2  xz     z^2    yz     xy\n");
        else if (maxorb == 3)
            std::printf("    Atom   Total    s     px     py     pz\n");
        else
            std::printf("    Atom   Total    s\n");
        std::printf("%8s Atomic\n", "");
        for (int i = 1; i <= numat; ++i) {
            double sum = 0.0;
            for (int j = nfirst[i]; j <= nlast[i]; ++j)
                sum += pa[(j * (j + 1)) / 2];
            std::string ssp = elemts_C::elemnt[nat[i]];
            while (!ssp.empty() && ssp.back() == ' ') ssp.pop_back();
            std::printf("%5d %2s%8.4f", i, ssp.c_str(), sum);
            for (int j = nfirst[i]; j <= nlast[i]; ++j)
                std::printf("%7.3f", pa[(j * (j + 1)) / 2]);
            std::printf("\n");
        }
        // F90 writmo.F90:805-815: HYPERFINE COUPLING COEFFICIENTS.
        // q(i) = pa_ii/3 + C_{i,nalpha-1}^2 * 2/3  (Fermi contact + spin density).
        if (keywrd.find(" HYPERFINE") != std::string::npos) {
            std::printf("\n\n%10s    HYPERFINE COUPLING COEFFICIENTS\n\n", "");
            int jj = nalpha - 1;
            if (jj < 1) jj = 1;
            std::vector<double> qhf((size_t)numat + 1, 0.0);
            for (int i = 1; i <= numat; ++i) {
                const int k0 = nfirst[i];
                qhf[i] = pa[(k0 * (k0 + 1)) / 2] * 0.3333333 +
                         c[k0][jj] * c[k0][jj] * 0.66666666;
            }
            std::string hline;
            for (int i = 1; i <= numat; ++i) {
                std::string sym = elemts_C::elemnt[nat[i]];
                while (!sym.empty() && sym.back() == ' ') sym.pop_back();
                char hb[64];
                std::snprintf(hb, sizeof(hb), "  %2s%2d%9.5f ", sym.c_str(), i, qhf[i]);
                hline += hb;
                if (i % 5 == 0) {
                    while (!hline.empty() && hline.back() == ' ') hline.pop_back();
                    std::printf("%s\n", hline.c_str());
                    hline.clear();
                }
            }
            if (!hline.empty()) {
                while (!hline.empty() && hline.back() == ' ') hline.pop_back();
                std::printf("%s\n", hline.c_str());
            }
        }
        for (size_t vi = 0; vi < p.size(); ++vi) pa[vi] = p[vi] - pb[vi];
    }
    // BONDS / ALLBO (F90 L819-832; molval M.O.-contribution part not ported).
    if (keywrd.find(" BONDS") != std::string::npos ||
        keywrd.find(" ALLBO") != std::string::npos) {
        bonds();
    }
    // 1ELE and ENPART (F90 L842-846).
    if (keywrd.find(" 1ELE") != std::string::npos) {
        std::printf(" FINAL ONE-ELECTRON MATRIX \n");
        vecprt(&h[0], norbs);
    }
    if (keywrd.find(" ENPART") != std::string::npos) enpart();
    // MULLIK / GRAPH (F90 L857-887; MOZYME lmo_to_eigenvectors part not ported).
    if ((keywrd.find(" MULLIK") != std::string::npos ||
         keywrd.find(" GRAPH") != std::string::npos) && !gui) {
        if (keywrd.find(" MULLIK") != std::string::npos)
            std::printf("\n%10s MULLIKEN POPULATION ANALYSIS\n", "");
        if (mozyme) {
            // Convert LMOs to eigenvectors (F90 writmo.F90 L857-887).
            if ((int)c.size() < norbs + 1)
                c.assign((size_t)norbs + 1, std::vector<double>((size_t)norbs + 1, 0.0));
            else
                for (int k = 1; k <= norbs; ++k)
                    if ((int)c[k].size() < norbs + 1) c[k].assign((size_t)norbs + 1, 0.0);
            lmo_to_eigenvectors(MOZYME_C::noccupied, MOZYME_C::ncf, MOZYME_C::nncf,
                                MOZYME_C::ncocc, MOZYME_C::noccupied,
                                MOZYME_C::icocc, MOZYME_C::icocc_dim,
                                MOZYME_C::cocc, MOZYME_C::cocc_dim, eigs, 0, c, 0);
            lmo_to_eigenvectors(MOZYME_C::nvirtual, MOZYME_C::nce, MOZYME_C::nnce,
                                MOZYME_C::ncvir, MOZYME_C::nvirtual,
                                MOZYME_C::icvir, MOZYME_C::icvir_dim,
                                MOZYME_C::cvir, MOZYME_C::cvir_dim, eigs,
                                MOZYME_C::noccupied, c, MOZYME_C::noccupied);
            // Convert "h" from MOZYME form to lower-half-triangle.
            const int mpack_l = (norbs * (norbs + 1)) / 2;
            std::vector<double> pb2((size_t)mpack_l + 1, 0.0);
            convert_mat_packed_to_triangle(h, pb2);
            h.assign((size_t)mpack_l + 1, 0.0);
            for (int k = 1; k <= mpack_l; ++k) h[k] = pb2[k];
        }
        mullik();
        if (keywrd.find(" GRAPH") != std::string::npos)
            std::printf("\n%10s DATA FOR GRAPH WRITTEN TO DISK\n", "");
    }
    // F90 writmo.F90 L890-894: SYBYL or Mulliken population printout.
    if (keywrd.find(" MULLIK") != std::string::npos)
        mpcpop(-1, nullptr);
    // ------------------------------------------------------------------
    // 专用输出块状态（2026-10-25 复核）
    // ------------------------------------------------------------------
    // 已移植/接线：MOZYME LMO 输出（prtlmo + eigs 能量）、MULLIK/GRAPH、
    // MOZYME 梯度表（prtgra）、周期体系 write_cell/write_pressure/
    // write_unit_cell_HOF、对称化 symtry（ndep!=0）、pKa（prtpka）。
    // 仍未接线（官方主 .out 一般不出现）：
    //   1) 超精细耦合常数（writmo.F90 ~L807，需超精细积分）；
    //   2) 公式/参数使用表（write_paras_used，官方仅 in_house_only 调用）；
    //   3) 少数 GUI/内部调试路径（modchg/modgra/web_message 等）。
    // 说明：PDB 输出由 run_mopac 的 PDBOUT 分支调用 pdbout（已接线）；
    // iwrite（archive 摘要文件）整段未移植，属独立大块。

    // ==================================================================
    // F90 writmo.F90:895-1148 — archive (.arc) summary file (iwrite = iarc).
    // F90 在 writmo 内部打开 archive_fn；C++ 用 fopen 实现，并临时把
    // chanel_C::out_fp 指向该文件，让 geout/wrttxt/dimens/timout 等既有
    // 格式化函数按 .arc 风格写同一份内容。
    // ==================================================================
    if (!gui) {
        std::FILE* afp = std::fopen(chanel_C::archive_fn.c_str(), "w");
        if (afp != nullptr) {
            chanel_C::out_fp = afp;
            // --- 摘要头（F90 L921-937） ---
            if (site_no > 9999)
                std::fprintf(afp, "\n\n%20s SUMMARY OF %s CALCULATION, Site No:%6d\n\n",
                             "", caltyp.c_str(), site_no);
            else if (site_no > -1)
                std::fprintf(afp, "\n\n%20s SUMMARY OF %s CALCULATION, Site No:%5d\n\n",
                             "", caltyp.c_str(), site_no);
            else
                std::fprintf(afp, "\n\n%20s SUMMARY OF %s CALCULATION\n\n",
                             "", caltyp.c_str());
            std::fprintf(afp, "%55sMOPAC2016 (Version: %s)\n", "", verson.c_str());
            std::fprintf(afp, "%55s%s\n", "", idate.c_str());
            if (!formula.empty()) std::fprintf(afp, "\n%s\n\n", formula.c_str());
            wrttxt(chanel_C::iarc);
            std::fprintf(afp, "\n\n\n    %-58.58s\n", flepo[iflepo]);
            std::fprintf(afp, "    %-58.58s\n", iter_msg[iscf]);
            // --- 能量块（F90 L941-987） ---
            if (stress > -0.5) {
                std::fprintf(afp, "\n\n\n%10sFINAL H.O.F PLUS STRESS =%17.5f KCAL/MOL =%14.5f KJ/MOL\n",
                             "", escf, escf * 4.184);
                std::fprintf(afp, "%10sFINAL STRESS            =%17.5f KCAL/MOL =%14.5f KJ/MOL\n",
                             "", stress, stress * 4.184);
                std::fprintf(afp, "%10sHEAT OF FORMATION       =%17.5f KCAL/MOL =%14.5f KJ/MOL\n",
                             "", escf - stress, (escf - stress) * 4.184);
            } else {
                std::fprintf(afp, "\n%10sHEAT OF FORMATION       =%17.5f KCAL/MOL =%14.5f KJ/MOL\n",
                             "", escf, escf * 4.184);
            }
            if (keywrd.find(" DISP") != std::string::npos) {
                std::fprintf(afp, "\n%10sTOTAL ENERGY            =%17.5f KCAL/MOL\n", "",
                             (elect + enuclr + solv_energy) * funcon_C::fpc_9);
                std::fprintf(afp, "%10sENERGY OF ATOMS         =%17.5f KCAL/MOL\n", "", atheat);
                std::fprintf(afp, "%10s                    SUM =%17.5f KCAL/MOL\n", "",
                             (elect + enuclr) * funcon_C::fpc_9 + atheat +
                                 solv_energy * funcon_C::fpc_9);
                if (std::fabs(hpress) > 1.0e-5)
                    std::fprintf(afp, "%10sENERGY DUE TO PRESSURE  =%17.5f KCAL/MOL\n", "", hpress);
                std::fprintf(afp, "%10sDISPERSION ENERGY       =%17.5f KCAL/MOL\n", "", E_disp);
                if (E_hb < -1.0e-5)
                    std::fprintf(afp, "%10sH-BOND ENERGY           =%17.5f KCAL/MOL\n", "", E_hb);
                if (std::fabs(nsp2_corr) > 1.0e-5)
                    std::fprintf(afp, "%10sMM CORRECTION FOR >N-   =%17.5f KCAL/MOL\n", "", nsp2_corr);
                if (std::fabs(Si_O_H_corr) > 1.0e-5)
                    std::fprintf(afp, "%10sMM CORR. FOR Si-O-H     =%17.5f KCAL/MOL\n", "", Si_O_H_corr);
                if (std::fabs(sum_dihed) > 1.0e-5)
                    std::fprintf(afp, "%10sMM CORR. FOR -CO-NH-    =%17.5f KCAL/MOL\n", "", sum_dihed);
                double ssum = (elect + enuclr) * funcon_C::fpc_9 + atheat + hpress +
                              solv_energy * funcon_C::fpc_9 + nsp2_corr + Si_O_H_corr +
                              sum_dihed + E_disp + E_hb;
                std::fprintf(afp, "%30sSUM =%17.5f KCAL/MOL\n\n", "", ssum);
                if (std::fabs(ssum - escf) > 1.0e-2)
                    std::fprintf(afp, " WARNING - An energy term is missing!\n");
                if (N_Hbonds > 0)
                    std::fprintf(afp, "%10sNo. OF HYDROGEN BONDS   =%11d       (H-bond Energy < -1.0 Kcal/mol)\n",
                                 "", N_Hbonds);
                if (E_hh > 1.0e-5)
                    std::fprintf(afp, "%10sH - H CORRECTION ENERGY =%17.5f KCAL/MOL\n", "", E_hh);
            }
            if (id == 3) write_unit_cell_HOF(chanel_C::iarc);
            const bool arch_state =
                (uhf || nopen != nclose || fract < 1.99) &&
                symmetry_C::state_Irred_Rep.find_first_not_of(' ') != std::string::npos;
            if (arch_state) {
                std::fprintf(afp, "%10sTOTAL ENERGY            =%17.5f EV   STATE:  %2d %s %s\n",
                             "", elect + enuclr + solv_energy, symmetry_C::state_QN,
                             symmetry_C::state_spin.c_str(),
                             symmetry_C::state_Irred_Rep.c_str());
            } else {
                std::fprintf(afp, "%10sTOTAL ENERGY            =%17.5f EV\n", "",
                             elect + enuclr + solv_energy);
            }
            std::fprintf(afp, "%10sELECTRONIC ENERGY       =%17.5f EV\n", "", elect);
            std::fprintf(afp, "%10sCORE-CORE REPULSION     =%17.5f EV\n", "", enuclr);
            if (cosmo_C::iseps)
                std::fprintf(afp, "%10sDIELECTRIC ENERGY       =%17.5f EV\n", "", cosmo_C::ediel);
            if (std::fabs(pressure) > 1.0e-4 && id == 3) {
                std::fprintf(afp, "%10sENERGY DUE TO PRESSURE  =%17.5f KCAL/MOL\n", "", hpress);
                const double vol3 = 0.0;  // 占位：周期体积已在主输出块给出
                (void)vol3;
            }
            if (lprtgra || gnorm > 1.0e-15)
                std::fprintf(afp, "%10sGRADIENT NORM           =%17.5f\n", "", gnorm);
            if (id == 0) {
                std::vector<double> dipv(4, 0.0);
                double dipv0 = 0.0;
                if (mozyme) dipv0 = dipole_for_MOZYME(dipv, 0);
                else dipv0 = dipole(p, coord, dipv, 0);
                if (mozyme)
                    std::fprintf(afp, "%10sDIPOLE                  =%17.5f DEBYE   \n",
                                 "", dipv0);
                else
                    std::fprintf(afp, "%10sDIPOLE                  =%17.5f DEBYE    POINT GROUP:       %-4s\n",
                                 "", dipv0, symmetry_C::name.c_str());
            }
            if (uhf) {
                std::fprintf(afp, "%10s(SZ)                    =  %16.6f\n", "", sz);
                if (fract < 1.0e-5 || fract > 0.9999)
                    std::fprintf(afp, "%10s(S**2)                  =  %16.6f\n", "", ss2);
                else
                    std::fprintf(afp, "%10s Average over configurations used, so (S**2) is not meaningful\n", "");
                std::fprintf(afp, "%10sNO. OF ALPHA ELECTRONS  =  %9d\n", "",
                             nalpha + (int)std::lround((nalpha_open - nalpha) * fract));
                std::fprintf(afp, "%10sNO. OF BETA  ELECTRONS  =  %9d\n", "",
                             nbeta + (int)std::lround((nbeta_open - nbeta) * fract));
            } else {
                const int iop_arc = (fract > 1.99) ? nopen : nclose;
                const int nopn_arc = nopen - iop_arc;
                std::fprintf(afp, "%10sNO. OF FILLED LEVELS    =  %9d\n", "",
                             nopen - nopn_arc);
                if (nopn_arc != 0)
                    std::fprintf(afp, "%10sAND NO. OF OPEN LEVELS  =  %9d\n", "",
                                 nopn_arc);
            }
            if (kchrge != 0)
                std::fprintf(afp, "%10sCHARGE ON SYSTEM        =  %9d\n", "", kchrge);
            if (!mozyme) {
                std::fprintf(afp, "%10sIONIZATION POTENTIAL    =  %16.6f EV\n", "", eionis);
                if (uhf) {
                    if (nalpha >= 1)
                        std::fprintf(afp, "%10sALPHA SOMO LUMO (EV)    =  %13.3f%7.3f\n",
                                     "", eigs[nalpha],
                                     (nalpha + 1 <= norbs) ? eigs[nalpha + 1] : 0.0);
                    if (nbeta >= 1)
                        std::fprintf(afp, "%10sBETA  SOMO LUMO (EV)    =  %13.3f%7.3f\n",
                                     "", eigb[nbeta],
                                     (nbeta + 1 <= norbs) ? eigb[nbeta + 1] : 0.0);
                } else if (nopen == nclose) {
                    if (nopen >= 1)
                        std::fprintf(afp, "%10sHOMO LUMO ENERGIES (EV) =  %13.3f%7.3f\n",
                                     "", eigs[nopen],
                                     (nopen + 1 <= norbs) ? eigs[nopen + 1] : 0.0);
                }
            }
            std::fprintf(afp, "%10sMOLECULAR WEIGHT        =  %14.4f\n", "", mol_weight);
            if (id == 3) write_cell(chanel_C::iarc);
            if (id == 3) write_pressure(chanel_C::iarc);
            if (cosmo_C::fepsi > 1.0e-10 && cosmo_C::area > 1.0e-3) {
                std::fprintf(afp, "%10sCOSMO AREA              =%14.2f SQUARE ANGSTROMS\n",
                             "", cosmo_C::area);
                std::fprintf(afp, "%10sCOSMO VOLUME            =%14.2f CUBIC ANGSTROMS\n",
                             "", cosmo_C::cosvol);
            }
            dimens(coord, chanel_C::iarc);
            std::fprintf(afp, "%10sSCF CALCULATIONS        =  %9d\n", "", nscf);
            timout(chanel_C::iarc);
            if (keywrd.find(" PRTCHAR") != std::string::npos)
                std::fprintf(afp, "\n\n%10sFINAL GEOMETRY OBTAINED%36s CHARGE\n", "", "");
            else
                std::fprintf(afp, "\n\n%10sFINAL GEOMETRY OBTAINED\n", "");
            geout(chanel_C::iarc);
            if (keywrd.find(" AIGOUT") != std::string::npos) {
                std::fprintf(afp, "\n\n  GEOMETRY IN GAUSSIAN Z-MATRIX STYLE\n");
                wrttxt(chanel_C::iarc);
                geoutg(chanel_C::iarc);
            }
            if (in_house_only) write_paras_used();
            std::fflush(afp);
            chanel_C::out_fp = nullptr;
            std::fclose(afp);
        }
    }
    ++numcal;
}
