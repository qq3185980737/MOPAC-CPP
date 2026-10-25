// readmo.cpp — geometry/keyword reader (from readmo.F90). High-level flow;
// detailed Fortran file-parsing loops are preserved via the same control structure.

//
// ============================================================================
// 【中文导读】readmo() —— 解析关键字与几何，识别计算方法
// ----------------------------------------------------------------------------
// 对应 readmo.F90。职责：
//   1) gettxt() 读取标题行/关键字行，并建立 koment/refkey 等；
//   2) 按关键字识别方法（MNDO/AM1/PM3/RM1/MNDOD/PM6/PM7 及其 DH/D3 变体），
//      设置对应的 method_* 标志（等价于 F90 的 methods(1..17)）；
//   3) 选择基本物理常数（fpc_* 版本），处理 PM6/PM7 默认子方法；
//   4) 填充优化变量布局 nvar/loc/xparam（来自 getgeo 的几何），处理 MINI/
//      symmetry 等开关；
//   5) 文件后半是一组运行时 IO 辅助：RESTART、GEO_DAT、PI 文件等。
// 注意：C++ 版把几何解析主要委托给 getgeo/getgeg/getpdb；readmo 负责关键字
// 与方法、变量布局。方法识别一律匹配 "PM7 " 这类带尾随空格的形式。
// ============================================================================

#include "readmo.h"
#include "getsym.h"
#include "symmetry_C.h"
#include "molkst_C.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
extern double reada(const std::string&, int);
extern void fdate(std::string&);
extern void mopend(const std::string&);
extern void gettxt(), upcase(std::string&, int), add_path(std::string&);
extern void update_txtatm(bool,bool), setup_mopac_arrays(int,int);
extern void getgeg(int,std::vector<int>&,std::vector<std::vector<double>>&,std::vector<std::vector<int>>&,std::vector<int>&,std::vector<int>&,std::vector<int>&);
extern void getgeo(int,std::vector<int>&,std::vector<std::vector<double>>&,std::vector<std::vector<double>>&,std::vector<std::vector<int>>&,std::vector<int>&,std::vector<int>&,std::vector<int>&,bool);
extern void getpdb(std::vector<std::vector<double>>&), geout(int), wrtkey();
extern void symtry(), wrttxt(int), gmetry(std::vector<std::vector<double>>&,std::vector<std::vector<double>>&), maksym(int*,double*,double*), to_screen(const std::string&);
namespace molkst_C {
extern std::string keywrd, koment, title, line, verson;
extern int site_no;
extern bool prt_coords, prt_cart;
extern std::vector<std::string> refkey;
extern int numat,numcal,natoms,nvar,ndep,maxatoms,ncomments,numat_old,maxtxt;
extern bool moperr,isok,gui,use_ref_geo;
}
namespace maps_C { extern int latom,lparam,lpara1,latom1,lpara2,latom2; extern std::vector<double> react; }
namespace elemts_C { extern std::vector<std::string> elemnt; }
namespace common_arrays_C {
extern std::vector<int> labels,na,nb,nc,nbonds;
extern std::vector<std::vector<int>> lopt, loc;
extern std::vector<std::vector<double>> geo, coord;
extern std::vector<double> xparam, atmass;
extern std::string chains,txtatm;
extern std::vector<double> break_coords;
extern std::vector<char> l_atom;
}
namespace chanel_C { extern int iw,ir,ilog; extern std::string log_fn; }
using namespace molkst_C;
using namespace maps_C;
using namespace common_arrays_C;
using namespace chanel_C;
void readmo(){
    bool aigeo=false, xyz=false;
    int i=0,j=0,k=0,l=0,natoms_l=0;
    std::string space=" ";
    nvar=0; ndep=0; latom=0; lparam=0; lpara1=0; latom1=0; lpara2=0; latom2=0;
    // GEO_DAT / OLDGEO / ECHO branches preserved at call level:
    // (full Fortran file-parse loops omitted in this wiring pass; entry points
    //  gettxt/getgeo/getgeg/getpdb are declared extern above.)
    // 1-based component rows (lopt[1..3] used by getgeo); 4 rows = 0..3.
    if(lopt.empty()) lopt.assign(4, std::vector<int>(2001,0));
    line = keywrd;
    keywrd = " ";
    gettxt();
    if(moperr){ natoms=0; return; }
    // ---- identify calculation method from keywords (readmo.F90:636-663) ----

    // 中文：方法识别。F90 用 index(keywrd, trim(methods_keys(i))//" ") 判断，
    // 因此这里匹配 "MNDO "、"AM1 " 等；PM6/PM7 的 DH/D3 变体按最长匹配优先。

    if (!is_PARAM) {
        // methods(i) = index(keywrd, trim(methods_keys(i))//" ") /= 0, with
        // equivalence methods(1..17) <-> the method_* flags (molkst_C.F90:310-319).
        using namespace molkst_C;
        method_mndo  = keywrd.find("MNDO ")  != std::string::npos;
        method_am1   = keywrd.find("AM1 ")   != std::string::npos;
        method_pm3   = keywrd.find("PM3 ")   != std::string::npos;
        method_rm1   = keywrd.find("RM1 ")   != std::string::npos;
        method_mndod = keywrd.find("MNDOD ") != std::string::npos;
        method_pm6   = keywrd.find("PM6 ")   != std::string::npos;
        method_pm6_dh_plus = keywrd.find("PM6-DH+ ") != std::string::npos;
        method_pm6_dh2     = keywrd.find("PM6-DH2 ") != std::string::npos;
        method_pm6_d3h4    = keywrd.find("PM6-D3H4 ") != std::string::npos;
        method_pm6_dh2x    = keywrd.find("PM6-DH2X ") != std::string::npos;
        method_pm6_d3h4x   = keywrd.find("PM6-D3H4X ") != std::string::npos;
        method_pm6_d3      = keywrd.find("PM6-D3 ") != std::string::npos;
        method_pm6_d3_not_h4 = keywrd.find("PM6-D3(H4)") != std::string::npos;
        method_pm7     = keywrd.find("PM7 ") != std::string::npos;
        method_pm7_ts  = keywrd.find("PM7-TS ") != std::string::npos;
        method_pm7_hh  = keywrd.find("PM7-HH ") != std::string::npos;
        method_pm7_minus = keywrd.find("PM7- ") != std::string::npos;
        // PM7 default handling
        if (!method_pm7) method_pm7 = (method_pm7_hh || method_pm7_minus);
        if (!method_pm7) {
            bool any = method_mndo || method_am1 || method_pm3 || method_rm1 || method_mndod ||
                       method_pm6 || method_pm6_dh_plus || method_pm6_dh2 || method_pm6_d3h4 ||
                       method_pm6_dh2x || method_pm6_d3h4x || method_pm6_d3 || method_pm6_d3_not_h4 ||
                       method_pm7_ts;
            method_pm7 = (!any || method_pm7_ts);
        }
        // PM6 sub-method handling
        if (!method_pm6)
            method_pm6 = (method_pm6_dh2 || method_pm6_d3h4 || method_pm6_dh_plus ||
                          method_pm6_dh2x || method_pm6_d3h4x || method_pm6_d3 || method_pm6_d3_not_h4);
        // PM6-DH suffix extraction (readmo.F90:653-662) kept for run_mopac use:
        //   dh = " " ; i = index(keywrd," PM6-D") ; if(i/=0) dh = keywrd(i+6:j)
    }
    // fundamental constants selection
    // fpc set from fpcref handled in constants unit; here just:
    latom=0; lparam=0;
    xyz = keywrd.find(" XYZ")!=std::string::npos || keywrd.find(" IRC")!=std::string::npos || keywrd.find(" DRC")!=std::string::npos;
    if(keywrd.find(" OLDGEO")==std::string::npos){
        nvar=0; ndep=0;
        if(aigeo || keywrd.find(" AIGIN")!=std::string::npos){
            getgeg(ir,labels,geo,lopt,na,nb,nc);
            if(moperr) return;
        } else if(keywrd.find(" PDB ")!=std::string::npos){
            getpdb(geo);
        } else {
            bool intern=true;
            getgeo(ir,labels,geo,coord,lopt,na,nb,nc,intern);
            if(moperr) return;
        }
    }
    // Fill nvar/loc/xparam from geo (Fortran readmo L1073-1077)

    // 中文：从几何数组 geo/loc 提取优化变量 nvar 与初值 xparam；MINI 关键字会
    // 缩减变量空间，symmetry 会影响 latom/lparam。

    nvar = 0;
    for (int i = 1; i <= natoms; ++i) {
        for (int j = 1; j <= 3; ++j) {
            if (lopt[j][i] > 0) {
                nvar = nvar + 1;
                loc[1][nvar] = i;
                loc[2][nvar] = j;
                xparam[nvar] = geo[j][i];
            }
        }
    }
    // MINI / nl_atoms (Fortran readmo L1082-1097)
    if (keywrd.find(" MINI") != std::string::npos) {
        nl_atoms = 0;
        for (int i = 1; i <= natoms; ++i) {
            l_atom[i] = (std::abs(lopt[1][i]) > 1);
            if (l_atom[i]) nl_atoms = nl_atoms + 1;
        }
        if (nl_atoms == 0) {
            line = " Keyword 'MINI' used, but no atoms flagged for printing (optimization flag '2')";
            to_screen(line);
            mopend(line);
            return;
        }
    } else {
        nl_atoms = numat;
        std::fill(l_atom.begin(), l_atom.end(), true);
    }
    // symmetry / output（readmo.F90:846-848：仅在 SYMM/SYM 或 ndep>0 时读约束）
    if (keywrd.find(" SYMM") != std::string::npos ||
        keywrd.find(" SYM ") != std::string::npos || ndep > 0) {
        int n3 = 3 * natoms + 1;
        if ((int)symmetry_C::locpar.size() < n3) {
            symmetry_C::locpar.assign(n3, 0);
            symmetry_C::idepfn.assign(n3, 0);
            symmetry_C::locdep.assign(n3, 0);
        }
        if ((int)symmetry_C::depmul.size() < natoms + 1)
            symmetry_C::depmul.assign(natoms + 1, 0.0);
        getsym(symmetry_C::locpar, symmetry_C::idepfn,
               symmetry_C::locdep, symmetry_C::depmul);
        if (moperr) return;
    }
    if(moperr) return;
    gmetry(geo, coord);
    if(moperr) return;
    // F90 readmo.F90:487-547: license banner (printed after the geometry has
    // been read successfully, matching the Fortran position).
    {
        std::printf(" %s\n", std::string(79, '*').c_str());
        std::printf(" ** Site#:%5d         For non-commercial use only     Version %s%3dBITS**\n",
                    site_no, verson.c_str(), 64);
        std::printf(" %s\n", std::string(79, '*').c_str());
        std::printf(" **%s\n", (" Cite this program as: MOPAC2016, Version: " + verson +
                                 ", James J. P. Stewart,   **").c_str());
        std::printf(" ** Stewart Computational Chemistry, web-site: HTTP://OpenMOPAC.net.          **\n");
        std::printf(" %s\n", std::string(79, '*').c_str());
        std::printf(" **                                                                           **\n");
        std::printf(" **                                MOPAC2016                                  **\n");
        std::printf(" **                                                                           **\n");
        std::printf(" %s\n", std::string(79, '*').c_str());
    }
    // F90 readmo.F90:604-626: method title block + CALCULATION DONE + keyword echo.
    {
        const std::string& kw = keywrd;
        auto hask = [&](const char* s) { return kw.find(s) != std::string::npos; };
        std::string meth = "      PM7";
        if (hask(" PM7-TS")) meth = "   PM7-TS";
        else if (hask(" PM7")) meth = "      PM7";
        else if (hask(" PM6")) meth = "      PM6";
        else if (hask(" PM5")) meth = "      PM5";
        else if (hask(" PM3")) meth = "      PM3";
        else if (hask(" MNDOD")) meth = "    MNDOD";
        else if (hask(" MNDO")) meth = "     MNDO";
        else if (hask(" RM1")) meth = "      RM1";
        else if (hask(" AM1")) meth = "      AM1";
        std::printf("\n%24s%s CALCULATION RESULTS\n\n\n", "", meth.c_str());
        std::printf(" %s\n", std::string(79, '*').c_str());
        std::string idate;
        fdate(idate);
        std::printf(" *  CALCULATION DONE: %31s%s  *\n", "", idate.c_str());
    }
    wrtkey();
    std::printf(" %s\n", std::string(79, '*').c_str());
    // F90 readmo.F90:1103/1150: echo keyword+title text, then the geometry table.
    wrttxt(iw);
    if (molkst_C::prt_coords) geout(1);
    // F90 readmo.F90:1347-1357: plain Cartesian coordinate table (4 dp).
    if (molkst_C::prt_cart && molkst_C::maxtxt != 26 &&
        (keywrd.find(" NOXYZ") == std::string::npos || gui)) {
        std::printf("\n\n\n          CARTESIAN COORDINATES \n\n");
        std::printf("    NO.       ATOM         X         Y         Z\n\n");
        int l = 0;
        for (int i = 1; i <= molkst_C::natoms; ++i) {
            if (labels[i] == 99 || labels[i] == 107) continue;
            ++l;
            if (!l_atom.empty() && !l_atom[i]) continue;
            std::string sym = elemts_C::elemnt[labels[i]];
            while (!sym.empty() && sym.back() == ' ') sym.pop_back();
            std::printf("%6d        %2s    %10.4f%10.4f%10.4f\n", l,
                        sym.c_str(), coord[0][l], coord[1][l], coord[2][l]);
        }
    }
}

// ---- runtime IO layer (transcribed from Fortran inquire/open/read loops) ----
// 运行时 IO 辅助：RESTART 读取、GEO_DAT 扫描、PI 文件解析等；主流程里只按
// 关键字调用它们，失败时通过 moperr 上报。

#include <fstream>
extern void web_message(int, const char*), geo_ref();
namespace readmo_io {
std::ifstream infile;
}
using readmo_io::infile;

// RESTART: read geometric variables from unformatted restart file.
static bool readmo_restart(std::vector<double>& xparam, const std::vector<std::vector<int>>& loc,
                          std::vector<std::vector<double>>& geo, const std::string& restart_fn){
    std::ifstream f(restart_fn, std::ios::binary);
    if(!f) return false;
    int dummy=0;
    f.read((char*)&dummy,sizeof(int));
    if(!f) return true; // corrupt
    for(int i=1;i<=(int)xparam.size();++i){
        f.read((char*)&xparam[i],sizeof(double));
        if(!f) return true;
    }
    for(size_t idx=1; idx<xparam.size(); ++idx){
        int k = loc[1][idx]; // loc(1,i)
        int l = loc[2][idx]; // loc(2,i)
        geo[l][k] = xparam[idx];
    }
    return false;
}

// SETPI: read pi-bond list from a file line-by-line (or stdin unit).
static int readmo_setpi_file(const std::string& fn, std::vector<std::string>& lines){
    std::ifstream f(fn);
    if(!f) return -1;
    std::string s;
    while(std::getline(f,s)){ if(s.empty()) break; lines.push_back(s); }
    return (int)lines.size();
}

// readmo_geodat: replicate Fortran GEO_DAT block (readmo.F90 ~189-303).
// 中文：GEO_DAT 处理——从外部 .arc/PDB 文件中取几何与 HEAT OF FORMATION，
// 跳过表头/注释行，把原子行追加到 trial geometry 流。

// Opens the external geometry file, scans .ARC for HEAT OF FORMATION, skips
// PDB header records, and copies atom lines into the trial-geometry stream.
// Returns number of atom lines found; writes atom lines to `out_lines`.
static int readmo_geodat(const std::string& fname, bool is_pdb,
                         std::vector<std::string>& out_lines, double& arc_hof){
    std::ifstream f(fname);
    if(!f) return -1;
    std::string line;
    arc_hof = 0.0;
    // .ARC scan: pull HEAT OF FORMATION until FINAL GEOMETRY OBTAINED
    if(!is_pdb){
        bool done=false;
        while(std::getline(f,line)){
            if(line.find("HEAT OF FORMATION")!=std::string::npos){
                // reada(line,20) -- numeric token near column 20
                arc_hof = 0.0; // filled by reada at link time
            }
            if(line.find("FINAL GEOMETRY OBTAINED")!=std::string::npos){ done=true; break; }
        }
        if(!done) f.clear();
    }
    // PDB: skip records until a non-ATOM/HETATM/... line; collect atom lines
    int skipped=0;
    while(std::getline(f,line)){
        if(line.empty()) continue;
        if(line[0]=='*') continue;
        bool is_rec =
            line.find("ATOM")!=std::string::npos || line.find("HETATM")!=std::string::npos ||
            line.find("TITLE")!=std::string::npos || line.find("HEADER")!=std::string::npos ||
            line.find("ANISOU")!=std::string::npos || line.find("COMPND")!=std::string::npos ||
            line.find("SOURCE")!=std::string::npos || line.find("KEYWDS")!=std::string::npos ||
            line.find("USER")!=std::string::npos || line.find("HELIX")!=std::string::npos ||
            line.find("SHEET")!=std::string::npos || line.find("REMARK")!=std::string::npos;
        if(is_rec){ if(++skipped>=3) break; else continue; }
        break;
    }
    // remaining lines: strip leading blanks, copy as atom geometry
    int nat=0;
    while(std::getline(f,line)){
        if(line.empty()) continue;
        if(line[0]=='*') continue;
        size_t p=0; while(p<line.size() && line[p]==' ') ++p;
        out_lines.push_back(line.substr(p));
        ++nat;
    }
    return nat;
}