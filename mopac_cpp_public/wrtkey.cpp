// wrtkey.cpp — C++ translation of "wrtkey.F90".
// Entry point: print recognized keyword groups and run conflict checks.
// wrtcon / wrtwor / wrtout are the three printer groups (large; stubbed).
// wrtchk: keyword-conflict validator (implemented from Fortran logic).

//
// ============================================================================
// 【中文导读】wrtkey() —— 关键字回显与冲突检查
// ----------------------------------------------------------------------------
// 对应 wrtkey.F90。打印本次作业使用的关键字分组（控制/几何/方法/输出），
// 并做关键字冲突校验（wrtchk）：例如多个方法互斥、C.I. 语法、MOZYME 家族、
// 非法组合等。
// 现状（重要）：wrtchk 已实现；三个打印分组 wrtcon/wrtwor/wrtout 目前是
// **空实现**（占位），因此 .out 中缺少关键字分组回显——属于已知缺口，
// 不影响计算，只影响输出完整性。
// ============================================================================

#include "wrtkey.h"

#include <algorithm>
#include <string>

#include "chanel_C.h"
#include "cosmo_C.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "mopend.h"
#include "reada.h"
#include "myword.h"
#include <cstdio>

using namespace chanel_C;
using namespace cosmo_C;
using namespace molkst_C;

// Printer groups (large keyword printers) — declared elsewhere / stubbed.

// wrtcon/wrtwor/wrtout：三个关键字打印分组——当前为空实现（只影响输出）。

struct WrtKeyEntry { const char* key; const char* text; };
static const WrtKeyEntry kCon[] = {
    {" MNDO ", " *  MNDO       - The MNDO Hamiltonian to be used"},
    {" AM1 ", " *  AM1        - The AM1 Hamiltonian to be used"},
    {" PM3 ", " *  PM3        - The PM3 Hamiltonian to be used"},
    {" MNDOD", " *  MNDOD      - The MNDOD Hamiltonian to be used"},
    {" PM6 ", " *  PM6        - The PM6 Hamiltonian to be used"},
    {" PM7-TS", " *  PM7-TS     - Calculate barrier height using PM7-TS"},
    {" PM7", " *  PM7        - The PM7 Hamiltonian to be used"},
    {" SPARKL", " *  SPARKLE    - Use SPARKLES when they exist."},
    {" RM1 ", " *  RM1        - The RM1 Hamiltonian to be used"},
    {" PM5 ", " *  PM5        - This Hamiltonian is not supported"},
    {"QMMM ", " *  QMMM       - Generate energies and gradients for use in MM codes"},
    {" BZ", " *  BZ         - Write file <name>.brz for use by program BZ"},
    {" BIRAD", " *  BIRADICAL- SYSTEM HAS TWO UNPAIRED ELECTRONS"},
    {" EXCI", " *  EXCITED    - FIRST EXCITED STATE IS TO BE OPTIMIZED"},
    {" VELO", " *  VELOCITY   - INPUT STARTING VELOCITIES FOR DRC"},
    {" GEO-OK", " *  GEO-OK     - OVERRIDE INTERATOMIC DISTANCE AND OTHER SAFETY CHECKS"},
    {" CHECK", " *  CHECK      - RUN EXTRA INTERATOMIC DISTANCE CHECKS"},
    {" AIGIN", " *  AIGIN      - GEOMETRY MUST BE IN GAUSSIAN FORMAT"},
    {" ESR", " *  ESR        - RHF SPIN DENSITY CALCULATION REQUESTED"},
    {" NOMM", " *  NOMM       - DO NOT MAKE MM CORRECTION TO CONH BARRIER"},
    {" MMOK", " *  MMOK       - APPLY MM CORRECTION TO CONH BARRIER"},
    {" CIS ", " *  CIS        - C.I. USES 1 ELECTRON EXCITATIONS ONLY"},
    {" CISD ", " *  CISD       - C.I. USES 1 AND 2 ELECTRON EXCITATIONS"},
    {" CISDT ", " *  CISDT      - C.I. USES 1, 2 AND 3 ELECTRON EXCITATIONS"},
    {" SING", " *  SINGLET    - SPIN STATE DEFINED AS A SINGLET"},
    {" DOUB", " *  DOUBLET    - SPIN STATE DEFINED AS A DOUBLET"},
    {" TRIP", " *  TRIPLET    - SPIN STATE DEFINED AS A TRIPLET"},
    {" QUAR", " *  QUARTET    - SPIN STATE DEFINED AS A QUARTET"},
    {" QUIN", " *  QUINTET    - SPIN STATE DEFINED AS A QUINTET"},
    {" SEXT", " *  SEXTET     - SPIN STATE DEFINED AS A SEXTET"},
    {" SEPT", " *  SEPTET     - SPIN STATE DEFINED AS A SEPTET"},
    {" OCTE", " *  OCTET      - SPIN STATE DEFINED AS A OCTET"},
    {" NONE", " *  NONET      - SPIN STATE DEFINED AS A NONET"},
    {" COSCCH", " *  COSCCH     - ADD IN COSMO CHARGE CORRECTIONS"},
    {" FIELD", " *  FIELD      - AN EXTERNAL ELECTRIC FIELD IS TO BE USED"},
    {" NOREOR", " *  NOREOR     - DO NOT REORIENTATE SYSTEM WHEN WORKING OUT THE SYMMETRY"},
    {" INVERT", " *  INVERT     - REVERSE ALL OPTIMIZATION FLAGS"},
    {" ESP ", " *  ESP        - ELECTROSTATIC POTENTIAL CALCULATION"},
    {" NSURF", " *  NSURF      - NUMBER OF LAYERS"},
    {" NOGPU", " *  NOGPU      - DO NOT USE GPU ACCELERATION"},
    {" SCALE", " *  SCALE      - SCALING FACTOR FOR VAN DER WAALS DISTANCE"},
    {" SCINCR", " *  SCINCR     - INCREMENT BETWEEN LAYERS"},
    {" SLOPE", " *  SLOPE      - SLOPE   - USED TO SCALE MNDO ESP CHARGES"},
    {" PDBOUT", " *  PDBOUT     - PRINT GEOMETY IN PDB FORMAT"},
    {" XENO", " *  XENO       - FRAGMENTS ARE ATTACHED TO RESIDUES"},
    {" CHARGES", " *  CHARGES    - IDENTIFY AND PRINT CHARGED ATOMS"},
    {" RESEQ", " *  RESEQ      - RESEQUENCE Z-MATIX INTO NORMAL PDB FORMAT"},
    {" NORES", " *  NORESEQ    - DO NOT RESEQUENCE Z-MATIX INTO NORMAL PDB FORMAT"},
    {" GEOCHK", " *  GEOCHK     - PRINT WORKING IN SUBROUTINE GEOCHK"},
    {" LEWIS", " *  LEWIS      - PRINT OUT LEWIS STRUCTURE, THEN STOP"},
    {" SETPI", " *  SETPI      - SOME OR ALL PI BONDS EXPLICITLY SET BY USER"},
    {" PDB ", " *  PDB        - INPUT GEOMETRY IS IN PDB FORMAT"},
    {" PDB(", " *  PDB(txt)   - SYMBOLS IN PDB FILE ARE DEFINED BY USER"},
    {" MOZ", " *  MOZYME     - USE LOCALIZED M.O.s IN SOLVING THE SCF EQUATIONS"},
    {" RAPID", " *  RAPID      - IN MOZYME SCF, USE ATOMS BEING OPTIMIZED ONLY"},
    {" REORTH", " *  REORTH     - RE-ORTHOGONALIZE LMOs EACH 10 SCF CALCULATIONS"},
    {" RE-LOCAL", " *  RE-LOCAL   - RE-LOCALIZE LMOs AT END OF CALCULATION"},
    {" SWAP", " *  SWAP       - THIS KEYWORD IS NOW REDUNDANT. SEE KEYWORD NOSWAP"},
    {" NOSWAP", " *  NOSWAP     - DO NOT SWAP ATOMS EVEN IF IT WILL IMPROVE OVERLAP IN GEO_REF"},
    {" A0 ", " *  A0         - INPUT GEOMETRY IS IN ATOMIC UNITS (A0)"},
    {" ANG", " *  ANGSTROMS- INPUT GEOMETRY IS IN ANGSTROMS"},
    {" ADD-H", " *  ADD-H      - ADD HYDROGEN ATOMS TO SATISFY VALENCE"},
    {" CONNOL", " *  CONNOLLY   - USE CONNOLLY SURFACE"},
    {" ESPRST", " *  ESPRST     - RESTART OF ELECTRIC POTENTIAL CALCULATION"},
    {" UHF", " *  UHF        - UNRESTRICTED HARTREE-FOCK CALCULATION"},
    {" RHF", " *  RHF        - RESTRICTED HARTREE-FOCK CALCULATION"},
    {" STATIC", " *  STATIC     - CALCULATE STATIC FIELD POLARIZABILITIES"},
    {" SETUP", " *  SETUP      - EXTRA KEYWORDS TO BE READ FROM FILE SETUP"},
    {" MAX", " *  MAX        - GRID SIZE 23*23 "},
    {" COSWRT", " *  COSWRT"},
    {" OLDCAV", " *  OLDCAV"},
    {" POLAR", " *  POLAR      - CALCULATE FIRST, SECOND AND THIRD-ORDER POLARIZABILITIES"},
    {" NORES", " *  NORES      - THIS IS THE DEFAULT.  USE RESIDUES IF RESIDUES ARE TO BE CALCULATED"},
};
static const WrtKeyEntry kOut[] = {
    {" PRTINT", " *  PRTINT     - INTERATOMIC DISTANCES TO BE PRINTED"},
    {" PRTCHAR", " *  PRTCHARGE- PRINT CHARGES IN ARC FILE"},
    {" MOPAC", " *  MOPAC      - USE OLD MOPAC CONVENTION FOR FIRST THREE ATOMS"},
    {" NOINT", " *  NOINTER    - THIS KEYWORD HAS BEEN REPLACED BY PRTINT"},
    {" ISOTOPE", " *  ISOTOPE    - FORCE MATRIX WRITTEN TO DISK (CHAN. 9 )"},
    {" DENOUT", " *  DENOUT     - DENSITY MATRIX OUTPUT ON CHANNEL 10"},
    {" OLDENS", " *  OLDENS     - INITIAL DENSITY MATRIX READ OF DISK"},
    {" CIOSCI", " *  CIOSCI     - PRINT WORKING IN SUBROUTINE CIOSCI"},
    {" ENPART", " *  ENPART     - ENERGY TO BE PARTITIONED INTO COMPONENTS"},
    {" NOXYZ", " *  NOXYZ      - CARTESIAN COORDINATES NOT TO BE PRINTED"},
    {" PRTXYZ", " *  PRTXYZ     - PRINT CARTESIAN COORDINATES"},
    {" NOTXT", " *  NOTXT      - DO NOT PRINT TEXT ASSOCIATED WITH AN ATOM"},
    {" GRAD", " *  GRADIENTS- ALL GRADIENTS TO BE PRINTED"},
    {" MINI", " *  MINI       - REDUCE OUTPUT BY ONLY PRINTING FLAGGED ATOMS"},
    {" VEC", " *  VECTORS    - FINAL MOLECULAR ORBITALS TO BE PRINTED"},
    {" VEC", " *  VECTORS    - FINAL EIGENVECTORS TO BE PRINTED"},
    {" CARTAB", " *  CARTAB     - PRINT ALL THE CHARACTER TABLES"},
    {" CHARST", " *  CHARST     - PRINT DETAILS OF SUBROUTINE CHARST"},
    {" DCART", " *  DCART      - PRINT DETAILS OF SUBROUTINE DCART"},
    {" DERI1", " *  DERI1      - PRINT DETAILS OF SUBROUTINE DERI1"},
    {" DERI2", " *  DERI2      - PRINT DETAILS OF SUBROUTINE DERI2 "},
    {" DERITR", " *  DERITR     - PRINT DETAILS OF SUBROUTINE DERITR"},
    {" DERIV", " *  DERIV      - PRINT DETAILS OF SUBROUTINE DERIV"},
    {" DERNVO", " *  DERNVO     - PRINT DETAILS OF SUBROUTINE DERNVO"},
    {" MAKVEC", " *  MAKVEC     - PRINT DETAILS OF SUBROUTINE MAKVEC"},
    {" DIIS", " *  DIIS       - PRINT DETAILS OF SUBROUTINE DIIS  "},
    {" FLEPO", " *  FLEPO      - PRINT DETAILS OF SUBROUTINE FLEPO "},
    {" FMAT", " *  FMAT       - PRINT DETAILS OF SUBROUTINE FMAT  "},
    {" DFORCE", " *  DFORCE     - PRINT FORCE MATRIX OVER CARTESIAN COORDINATES"},
    {" HCORE", " *  HCORE      - PRINT DETAILS OF SUBROUTINE HCORE "},
    {" MOLDAT", " *  MOLDAT     - PRINT DETAILS OF SUBROUTINE MOLDAT"},
    {" FREQCY", " *  FREQCY     - PRINT DETAILS OF SUBROUTINE FREQCY"},
    {" ITER", " *  ITER       - PRINT DETAILS OF SUBROUTINE ITER  "},
    {" LINMIN", " *  LINMIN     - PRINT DETAILS OF SUBROUTINE LINMIN"},
    {" MOLSYM", " *  MOLSYM     - PRINT DETAILS OF SUBROUTINE MOLSYM"},
    {" RAMA", " *  RAMA       - PRINT RAMACHANDRA ANGLES FOR PROTEIN RESIDUES"},
    {" SYMOIR", " *  SYMOIR     - PRINT DETAILS OF SUBROUTINE SYMOIR"},
    {" GROUP", " *  GROUP      - PRINT DETAILS OF SUBROUTINE GROUP "},
    {" SYMTRZ", " *  SYMTRZ     - PRINT DETAILS OF SUBROUTINE SYMTRZ"},
    {" DENS", " *  DENSITY    - FINAL DENSITY MATRIX TO BE PRINTED"},
    {" SPIN", " *  SPIN       - FINAL UHF SPIN MATRIX TO BE PRINTED"},
    {" PRNT", " *  PRNT       - EXTRA PRINTING IN EF ROUTINE"},
    {" DISP", " *  DISP       - PRINT DISPERSION AND HYDROGEN BOND ENERGIES"},
    {" TIMES", " *  TIMES      - TIMES OF VARIOUS STAGES TO BE PRINTED"},
    {" BONDS", " *  BONDS      - NON-HYDROGEN BOND-ORDERS TO BE PRINTED"},
    {" ALLBO", " *  ALLBONDS   - ALL SIGNIFICANT BOND-ORDERS TO BE PRINTED"},
    {" BONDS", " *  BONDS      - FINAL BOND-ORDER MATRIX TO BE PRINTED"},
    {" SYBYL", " *  SYBYL      - OUTPUT SYBYL FILE"},
    {" FOCK", " *  FOCK       - LAST FOCK MATRIX TO BE PRINTED"},
    {" LARGE", " *  LARGE      - EXPANDED OUTPUT TO BE PRINTED"},
    {" AIGOUT", " *  AIGOUT     - IN ARC FILE, INCLUDE AB-INITIO GEOMETRY"},
    {" GRAP", " *  GRAPH      - GENERATE FILE FOR GRAPHICS"},
    {" 1ELEC", " *  1ELECTRON- FINAL ONE-ELECTRON MATRIX TO BE PRINTED"},
    {" INTERP", " *  INTERP     - PRINT DETAILS OF CAMP-KING CONVERGER"},
    {" HTML", " *  HTML       - WRITE HTML SCRIPT TO READ PDB FILE USING JSMOL"},
    {" SUPER", " *  SUPER      - PRINT SUPERDELOCALIZABILITIES"},
    {" ALLVEC", " *  ALLVEC     - PRINT ALL EIGENVECTORS"},
    {" H-PRIO", " *  H-PRIOR    - HEAT OF FORMATION TAKES PRIORITY IN DRC"},
    {" X-PRIO", " *  X-PRIOR    - GEOMETRY CHANGES TAKE PRIORITY IN DRC"},
    {" T-PRIO", " *  T-PRIOR    - TIME TAKES PRIORITY IN DRC"},
    {" MECI", " *  MECI       - M.E.C.I. WORKING TO BE PRINTED"},
    {" HESSIAN", " *  HESSIAN    - WRITE OUT HESSIAN FROM GEOMERY OPTIMIZATION"},
    {" LOCAL", " *  LOCALIZE   - LOCALIZED ORBITALS TO BE PRINTED"},
    {" MULLIK", " *  MULLIK     - THE MULLIKEN ANALYSIS TO BE PERFORMED"},
    {" PI ", " *  PI         - BONDS MATRIX, SPLIT INTO SIGMA-PI-DELL COMPONENTS, TO BE PRINTED"},
    {" ECHO", " *  ECHO       - ALL INPUT DATA TO BE ECHOED BEFORE RUN"},
    {" DEBUG ", " *  DEBUG      - DEBUG OPTION TURNED ON"},
    {" POTWRT", " *  POTWRT     - WRITE OUT ELECTRIC POT. DATA TO FILE 21"},
    {" PMEP", " *  PMEP       - GENERATE THE WANG-FORD ELECTROSTATIC POTENTIAL"},
    {" PRTMEP", " *  PRTMEP     - MEP CONTOUR DATA OUTPUT TO <FILENAME>.mep"},
    {" MINMEP", " *  MINMEP     - PRINT MEP MINIMA IN THE PLANE DEFINED"},
    {" PL", " *  PL         - MONITOR CONVERGENCE IN DENSITY MATRIX"},
    {" PKA", " *  PKA        - PRINT pKa FOR MOST ACIDIC HYDROGEN"},
    {" POPS", " *  POPS       - PRINT SCF ATOMIC ORBITAL POPULATIONS"},
    {" BCC", " *  BCC        - THE SYSTEM IS BODY-CENTERED CUBIC"},
    {" EIGS", " *  EIGS       - PRINT ALL EIGENVALUES IN ITER"},
    {" HYPERF", " *  HYPERFINE- HYPERFINE COUPLING CONSTANTS TO BE PRINTED"},
    {" THERMO", " *  THERMO     - THERMODYNAMIC QUANTITIES TO BE CALCULATED"},
};
static const WrtKeyEntry kWor[] = {
    {" EIGINV", " *  EIGINV     - USE HESSIAN EIGENVALUE REVERSION IN EF"},
    {" NONR", " *  NONR       - DO NOT USE NEWTON-RAPHSON STEP IN EF"},
    {" SNAP", " *  SNAP       - INCREASE PRECISION OF SYMMETRY ANGLES"},
    {" CAMP", " *  CAMP,KING- THE CAMP-KING CONVERGER TO BE USED"},
    {" KING", " *  CAMP,KING- THE CAMP-KING CONVERGER TO BE USED"},
    {" LET", " *  LET        - OVERRIDE SOME SAFETY CHECKS"},
    {" OLDGEO", " *  OLDGEO     - PREVIOUS GEOMETRY TO BE USED"},
    {" OLDFPC", " *  OLDFPC     - OLD FUNDAMENTAL PHYSICAL CONSTANTS TO BE USED"},
    {" OLD_HESS", " *  OLD_HESS   - USE THE OLD HESSIAN MATRIX"},
    {" PREC", " *  PRECISE    - CRITERIA TO BE INCREASED BY 100 TIMES"},
    {" NOANCI", " *  NOANCI     - DO NOT USE ANALYTICAL C.I. DERIVATIVES"},
    {" DFP", " *  DFP        - USE DAVIDON FLETCHER POWELL OPTIMIZER"},
    {" XYZ", " *  XYZ        - CARTESIAN COORDINATE SYSTEM TO BE USED"},
    {" RESTART", " *  RESTART    - CALCULATION RESTARTED"},
    {" RSCAL", " *  RSCAL      - SCALE P-RFO STEP IN EF TO TRUST RADIUS"},
    {" SYMAVG", " *  SYMAVG     - AVERAGE SYMMETRY EQUIVALENT ESP CHARGES"},
    {" STO3G", " *  STO3G      - DEORTHOGONALIZE ORBITALS IN STO-3G BASIS"},
    {" FORCE ", " *  FORCE      - FORCE CALCULATION SPECIFIED"},
    {" FORCETS", " *  FORCETS    - VERIFY THAT TRANSITION STATE IS GENUINE"},
    {" BFGS", " *  BFGS       - USE THE BFGS GEOMETRY OPTIMIZER"},
    {" LBFGS", " *  LBFGS      - USE THE LBFGS GEOMETRY OPTIMIZER"},
    {" EF", " *  EF         - USE EF ROUTINE FOR MINIMUM SEARCH"},
    {" TS", " *  TS         - USE EF ROUTINE FOR TS SEARCH"},
    {" LOCATE-TS", " *  LOCATE-TS  - LOCATE A TRANSITION STATE IN A PROTEIN MECHANISM"},
    {" REFINE-TS", " *  REFINE-TS  - REFINE A TRANSITION STATE IN A PROTEIN MECHANISM"},
    {" NOSYM", " *  NOSYM      - POINT-GROUP SYMMETRY SET TO C1"},
    {" SMOOTH", " *  SMOOTH     - IN A GRID CALCULATION, REMOVE COMPUTATIONAL ARTIFACTS"},
    {" AUTOSYM", " *  AUTOSYM    - SYMMETRY TO BE IMPOSED AUTOMATICALLY"},
    {" 1SCF ", " *  1SCF       - DO 1 SCF AND THEN STOP "},
    {" SIGMA", " *  SIGMA      - GEOMETRY TO BE OPTIMIZED USING SIGMA."},
    {" NLLSQ", " *  NLLSQ      - GRADIENTS TO BE MINIMIZED USING NLLSQ."},
    {" SADDLE", " *  SADDLE     - TRANSITION STATE TO BE OPTIMIZED"},
    {" 0SCF", " *  0SCF       - AFTER READING AND PRINTING DATA, STOP"},
    {" QPMEP ", " *  QPMEP      - CHARGES DERIVED FROM WANG-FORD TYPE AM1 MEP"},
    {" NEWGEO", " *  NEWGEO     - USE BOTH INT AND CARTESIAN COORDINATES"},
    {" BIGSCF", " *  BIGSCF     - DO INITIAL FULL SCF WHEN RESTARTING JOB"},
    {" WILLIA", " *  WILLIAMS   - USE WILLIAMS SURFACE"},
    {" INT ", " *  INT        - INTERNAL COORDINATE SYSTEM TO BE USED"},
    {" PECI", " *  PECI       - SINGLE AND PAIRED ELECTRON EXCITATIONS ONLY"},
};

static void print_kw_group(const WrtKeyEntry* tab, size_t n, const std::string& allkey) {
    // F90 keywrd always carries a leading/trailing blank; C++ keywrd may not,
    // so pad to keep keys like " PM7" / " 1SCF " matchable.
    std::string ak = " " + allkey + " ";
    for (size_t i = 0; i < n; ++i)
        if (myword(ak, std::string(tab[i].key))) std::printf("%s\n", tab[i].text);
}

static void wrtcon(std::string& allkey) {
    allkey = keywrd;  // F90 wrtcon: allkey = trim(keywrd)
    print_kw_group(kCon, sizeof(kCon) / sizeof(kCon[0]), allkey);
}
static void wrtout(std::string& allkey) {
    print_kw_group(kOut, sizeof(kOut) / sizeof(kOut[0]), allkey);
}
static void wrtwor(std::string& allkey) {
    print_kw_group(kWor, sizeof(kWor) / sizeof(kWor[0]), allkey);
    // F90 wrtkey.F90:1173-1250: T= and DUMP=N time-limit lines.
    std::string ak = allkey;
    double tleft = 172800.0;  // two days
    std::string chrono = "SECONDS";
    double tmul = 1.0;
    if (myword(ak, " T=")) {
        int ip = (int)keywrd.find(" T=");
        tleft = reada(keywrd, ip + 1);
        std::string tok = keywrd.substr(ip + 3);
        size_t sp = tok.find(' ');
        if (sp != std::string::npos) tok = tok.substr(0, sp);
        if (tok.find('M') != std::string::npos) { chrono = "MINUTES"; tmul = 60.0; }
        else if (tok.find('H') != std::string::npos) { chrono = "HOURS"; tmul = 3600.0; }
        else if (tok.find('D') != std::string::npos) { chrono = "DAYS"; tmul = 86400.0; }
        else if (tok.find('W') != std::string::npos) { chrono = "WEEKS"; tmul = 7.0 * 86400.0; }
    }
    const bool had_t = keywrd.find(" T=") != std::string::npos;
    if (had_t && tleft >= 99999.9)
        std::printf(" *  T=         - A TIME OF%9.1g %s REQUESTED\n", tleft, chrono.c_str());
    else
        std::printf(" *  T=         - A TIME OF%9.1f %s REQUESTED\n", tleft, chrono.c_str());
    double tdump = 7200.0;
    chrono = "SECONDS";
    if (myword(ak, " DUMP")) {
        int ip = (int)keywrd.find(" DUMP");
        tdump = reada(keywrd, ip + 1);
        std::string tok = keywrd.substr(ip + 6);
        size_t sp = tok.find(' ');
        if (sp != std::string::npos) tok = tok.substr(0, sp);
        if (tok.find('M') != std::string::npos) chrono = "MINUTES";
        else if (tok.find('H') != std::string::npos) chrono = "HOURS";
        else if (tok.find('D') != std::string::npos) chrono = "DAYS";
        else if (tok.find('W') != std::string::npos) chrono = "WEEKS";
    }
    if (tdump < 99999.9)
        std::printf(" *  DUMP=N     - RESTART FILE WRITTEN EVERY%10.3f %s\n", tdump, chrono.c_str());
    else
        std::printf(" *  DUMP=N     - RESTART FILE WRITTEN EVERY%10.3g %s\n", tdump, chrono.c_str());
}


static bool has(const std::string& s, const std::string& w) {
    return s.find(w) != std::string::npos;
}

static void wrtchk(std::string& allkey) {
    (void)allkey;
    bool birad = has(keywrd, " BIRAD");
    bool exci = has(keywrd, " EXCI");
    bool ci = has(keywrd, " C.I.");
    bool trip = has(keywrd, " TRIP");
    uhf = has(keywrd, " UHF");
    rhf = has(keywrd, " RHF") || has(keywrd, " MECI") || ci ||
          (!uhf && has(keywrd, " OPEN"));

    // F90 wrtkey.F90:652-680: parse/validate the C.I. keyword.  A bare "C.I."
    // (no '=n' or '=(n1,n2)') is an error; otherwise set nmos exactly as the
    // Fortran does (nmos is consumed later by meci/densit).  Without this, an
    // "C.I.=n" job would run with nmos=0 and give a wrong CI energy.
    if (ci) {
        if (keywrd.find(" C.I.=(") != std::string::npos) {
            // C.I.=(n1,n2): nmos = Nint(reada(keywrd, Index("C.I.=(")+5))
            meci_C::nmos = (int)(reada(keywrd, (int)keywrd.find("C.I.=(") + 5) + 0.5);
        } else if (keywrd.find(" C.I.=") != std::string::npos) {
            // C.I.=n: nmos = Nint(reada(keywrd, Index("C.I.=")+5))
            meci_C::nmos = (int)(reada(keywrd, (int)keywrd.find("C.I.=") + 5) + 0.5);
        } else {
            mopend("C.I. keyword must be of form 'C.I.=n' or 'C.I.=(n1,n2)'");
            return;
        }
        if (meci_C::nmos > 29) {
            mopend("Maximum size of open space = 29 M.O.s");
            return;
        }
    }

    // F90 wrtkey.F90:303-306: MOZYME is selected when any MOZYME-family
    // keyword is present.  Set here (at keyword-echo time), independent of
    // compfg's iseps-guarded block, exactly as the Fortran does.
    {
        bool l_temp = (keywrd.find(" MOZ") != std::string::npos ||
                       keywrd.find(" LEWIS") != std::string::npos ||
                       keywrd.find(" LOCATE-TS") != std::string::npos ||
                       keywrd.find(" RESEQ") != std::string::npos ||
                       keywrd.find(" CHARGES") != std::string::npos ||
                       keywrd.find(" REFINE-TS") != std::string::npos ||
                       keywrd.find(" RAPID") != std::string::npos ||
                       keywrd.find(" SITE=") != std::string::npos ||
                       keywrd.find(" SITE(") != std::string::npos);
        mozyme = (l_temp || keywrd.find(" PDBOUT") != std::string::npos);
    }

    if (!(mozyme || has(keywrd, " PDBOUT") || has(keywrd, " RESID") ||
          has(keywrd, " ADD-H")) &&
        !has(keywrd, " 0SCF")) {
        if (has(keywrd, " PDBOUT")) {
            mopend("Keyword PDBOUT only works when MOZYME or 0SCF is also present");
            return;
        }
        if (has(keywrd, " RESID")) {
            mopend("Keyword RESIDUES only works when MOZYME or 0SCF is also present");
            return;
        }
        if (has(keywrd, " CVB")) {
            mopend("Keyword CVB only works with MOZYME");
            return;
        }
        if (has(keywrd, " SETPI")) {
            mopend("Keyword SETPI only works with MOZYME");
            return;
        }
    }
    if (mozyme) {
        if (uhf) mopend("Keyword UHF cannot be used with MOZYME");
        if (has(keywrd, " ENPART")) mopend("Keyword ENPART is not available with MOZYME");
        if (has(keywrd, " LOCAL")) mopend("Keyword LOCAL is not available with MOZYME");
        if (id != 0 && has(keywrd, " CUTOF")) {
            mopend("CUTOFx=n.nn -type keywords do not work with MOZYME for infinite systems");
            return;
        }
        if (has(keywrd, " 1SCF") && has(keywrd, " RAPID")) {
            mopend("RAPID cannot be used with 1SCF");
            return;
        }
    }
    if (has(keywrd, " MULLIK") && uhf) {
        mopend("MULLIKEN POPULATION NOT AVAILABLE WITH UHF");
        return;
    }
    if (uhf) {
        if (rhf) mopend("UHF and RHF cannot both be used");
        if (birad || exci || ci) {
            std::printf("UHF USED WITH EITHER BIRAD, EXCITED OR C.I.\n");
            mopend("IMPOSSIBLE OPTION REQUESTED");
            return;
        }
        if (has(keywrd, " POLAR")) {
            mopend("POLAR does not work with UHF");
            return;
        }
    } else if (exci && trip) {
        std::printf("EXCITED USED WITH TRIPLET\n");
        mopend("IMPOSSIBLE OPTION REQUESTED");
        return;
    }
    if (has(keywrd, " PMEP") && !method_am1) {
        mopend("PMEP only works with AM1");
        return;
    }
    if (has(keywrd, " INT ") && has(keywrd, " XYZ")) {
        mopend("INT cannot be used with XYZ");
        return;
    }
    if (id > 0 && has(keywrd, "EPS=")) {
        if (id == 1) std::printf("COSMO cannot be used with polymers\n");
        if (id == 2) std::printf("COSMO cannot be used with layer systems\n");
        if (id == 3) std::printf("COSMO cannot be used with solids\n");
        mopend("COSMO cannot be used with systems with Tv");
        return;
    }
    if (has(keywrd, " T-PRIO") && !has(keywrd, " DRC")) {
        std::printf("T-PRIO AND NO DRC\n");
        mopend("IMPOSSIBLE OPTION REQUESTED");
        return;
    }
    // Only one method allowed.
    int m = 0;
    if (method_am1) ++m;
    if (method_pm3) ++m;
    if (method_pm6) ++m;
    if (method_pm7) ++m;
    if (method_mndo) ++m;
    if (method_mndod) ++m;
    if (method_rm1) ++m;
    if (m > 1) {
        mopend("ONLY ONE OF MNDO, MNDOD, AM1, PM3, RM1, AND PM6 ALLOWED");
        return;
    }
    // Only one geometry option: strip quoted strings first.
    line = keywrd;
    for (size_t pos = 0; pos < line.size(); ++pos) {
        if (line[pos] == '"') {
            size_t e = line.find('"', pos + 1);
            if (e == std::string::npos) break;
            for (size_t k = pos; k <= e; ++k) line[k] = ' ';
        }
    }
    int g = 0;
    if (has(line, " BFGS")) ++g;
    if (has(line, " LBFGS")) ++g;
    if (has(line, " EF")) ++g;
    if (has(line, " TS")) ++g;
    if (has(line, " SIGMA")) ++g;
    if (has(line, " NLLSQ")) ++g;
    if (has(line, " FORCE") || has(line, " IRC") || has(line, " DRC")) ++g;
    if (g > 1) {
        mopend("MORE THAN ONE GEOMETRY OPTION HAS BEEN SPECIFIED. CONFLICT MUST BE RESOLVED BEFORE JOB WILL RUN.");
        return;
    }
    if (has(keywrd, " HESSIAN") && !has(keywrd, " EF")) {
        mopend("Keyword EF must be present if HESSIAN is used");
        return;
    }
}

// wrtkey(): 主入口——按关键字分组回显并调用 wrtchk 做冲突检查。

void wrtkey() {
    wrtcon(allkey);
    if (moperr) return;
    wrtwor(allkey);
    wrtout(allkey);
    // F90 wrtkey.F90 SOLVATION KEYWORDS block: assign fepsi from EPS= / PKA
    // (or the bare COSMO keyword -> EPS=999) and echo the EPS line.
    if (has(keywrd, " EPS=") || has(keywrd, " PKA")) {
        if (has(keywrd, "EPS=CRSDEF")) {
            fepsi = 1.0;
        } else {
            double epsi;
            if (has(keywrd, " PKA"))
                epsi = 78.4;
            else
                epsi = reada(keywrd, (int)keywrd.find("EPS="));
            std::printf(" *  EPS=%6.2f - USE ANDREAS KLAMT'S COSMO IMPLICIT SOLVATION MODEL\n", epsi);
            fepsi = (epsi - 1.0) / (epsi + 0.5);
            if (epsi < 1.0) {
                std::printf(" *  \n");
                std::printf(" *  The lower bound of a dielectric is 1.00 = vacuum\n");
                std::printf(" *  Keywords: %s\n", keywrd.c_str());
                std::printf(" *  \n");
                mopend("Nonsense value of dielectric constant supplied");
                return;
            }
        }
    } else if (has(keywrd, " COSMO")) {
        double epsi = 999.0;
        std::printf(" *  EPS=%6.2f - USE ANDREAS KLAMT'S COSMO IMPLICIT SOLVATION MODEL\n", epsi);
        fepsi = (epsi - 1.0) / (epsi + 0.5);
    }
    // F90 wrtkey.F90 752-758: NSPA default is 42 unless given explicitly.
    if (has(keywrd, " NSPA")) {
        size_t p = keywrd.find(" NSPA");
        int i = (int)(reada(keywrd, (int)p) + 0.5);
        std::printf(" *  NSPA=%5d\n", i);
        nspa = i;
    } else {
        nspa = 42;
    }
    // Fortran wrtkey L921-958: OUTPUT keyword controls the prt_* flags.
    // Without OUTPUT all print flags default to TRUE.
    if (has(keywrd, " OUTPUT")) {
        size_t i = keywrd.find(" OUTPUT");
        size_t j = keywrd.find(") ", i);
        if (j != std::string::npos && keywrd.find("OUTPUT(") != std::string::npos) {
            size_t a = keywrd.find("OUTPUT(") + 6;
            std::string line_ = keywrd.substr(a, j - a);
            prt_coords    = line_.find("C") != std::string::npos;
            prt_gradients = line_.find("G") != std::string::npos;
            prt_cart      = line_.find("X") != std::string::npos;
            prt_charges   = line_.find("Q") != std::string::npos;
            prt_pops      = line_.find("P") != std::string::npos;
            prt_topo      = line_.find("T") != std::string::npos;
        } else {
            prt_coords = false;
            prt_gradients = false;
            prt_cart = false;
            prt_charges = false;
            prt_pops = false;
            prt_topo = false;
        }
    } else {
        prt_coords = true;
        prt_gradients = true;
        prt_cart = true;
        prt_charges = true;
        prt_pops = true;
        prt_topo = true;
    }
    wrtchk(allkey);
}
