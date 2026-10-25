// ============================================================================
// 【中文导读】write_trajectory.cpp —— 输出/文件写出
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 write_trajectory.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// write_trajectory.cpp — C++ translation of "write_trajectory.F90".
// Append a geometry frame to the XYZ trajectory file (mode=1), or read back
// the XYZ file and rewrite the path in reverse (mode=2). Optional PDBOUT.

#include "write_trajectory.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"
#include "pdbout.h"
#include "to_screen.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;

void write_trajectory(double* xyz, int mode, double* charge, double escf,
                      double ekin, double time, double xtot) {
    static int icalcn = -1;
    static int imodel = 0;
    static int npt = 0;
    const int ipdb = 14;

    std::ofstream xyzout;
    std::ifstream xyzin;
    std::ofstream pdboutf;

    if (icalcn != step_num) {
        xyzout.open(xyz_fn, std::ios::out | std::ios::trunc);
        if (keywrd.find(" PDBOUT") != std::string::npos) {
            std::string pfn = xyz_fn.substr(0, xyz_fn.size() - 3) + "pdb";
            pdboutf.open(pfn);
        }
        icalcn = step_num;
    } else if (mode == 1) {
        xyzout.open(xyz_fn, std::ios::out | std::ios::app);
    }

    if (mode == 1) {
        // Write out a trajectory frame.
        xyzout << std::setw(6) << nl_atoms << " \n";
        int jloop = itemp_1;
        int nd1 = std::max(1, (int)std::log10(jloop * 1.01) + 1);
        double factor = std::fabs(escf);
        int nd2 = std::max(1, (int)std::log10(factor));
        int nd3 = std::max(1, (int)std::log10(4.184 * factor));
        {
            std::ostringstream os;
            os << "Profile." << std::setw(nd1) << jloop
               << " HEAT OF FORMATION =" << std::setw(nd2 + 9) << std::fixed
               << std::setprecision(3) << escf << " KCAL ="
               << std::setw(nd3 + 9) << std::fixed << std::setprecision(3)
               << escf * 4.184 << " KJ";
            line = os.str();
        }
        size_t kpos = line.find("HEAT");
        ++npt;
        xyzout << std::setw(8) << line.substr(0, 8) << std::setw(4) << npt
               << " " << line.substr(kpos) << "\n";
        for (int i = 1; i <= numat; ++i) {
            if (l_atom[i]) {
                // F90 format: (3X,A2,3F15.5)，且 Fortran 的 elemnt 是右对齐
                // 的 " H"（单字母元素），C++ 表是 "H "，这里按右对齐 2 字符输出。
                std::string sym = elemnt[nat[i]];
                size_t lastnz = sym.find_last_not_of(' ');
                if (lastnz != std::string::npos) sym.erase(lastnz + 1);
                else sym.clear();
                xyzout << "   " << std::setw(2) << sym << std::setw(15)
                       << std::fixed << std::setprecision(5)
                       << xyz[(i - 1) * 3 + 1] << std::setw(15)
                       << xyz[(i - 1) * 3 + 2] << std::setw(15)
                       << xyz[(i - 1) * 3 + 3] << "\n";
            }
        }
        if ((jloop % 10) == 0) {
            if (time > 1e-6) {
                std::ostringstream os;
                os << " CYCLE:" << std::setw(6) << jloop << "  Pot.E:"
                   << std::setw(17) << std::fixed << std::setprecision(5) << escf
                   << "  Kin.E:" << std::setw(9) << std::fixed
                   << std::setprecision(4) << ekin << "  Move:"
                   << std::setw(9) << std::fixed << std::setprecision(4) << xtot
                   << "  Time:" << std::setw(9) << std::fixed
                   << std::setprecision(3) << std::min(time, 99999.999);
                to_screen(os.str());
            } else {
                std::ostringstream os;
                os << " CYCLE:" << std::setw(6) << jloop
                   << "  Potential energy:" << std::setw(17) << std::fixed
                   << std::setprecision(5) << escf << "  Diff.:" << std::setw(9)
                   << std::fixed << std::setprecision(4) << ekin << "  Move:"
                   << std::setw(9) << std::fixed << std::setprecision(4) << xtot;
                to_screen(os.str());
            }
        }
        if (keywrd.find(" PDBOUT") != std::string::npos) {
            ++imodel;
            pdboutf << "MODEL " << std::setw(6) << imodel << "\n";
            pdbout(ipdb);
            pdboutf << "ENDMDL\n";
        }
        if (charge && charge[0] > -100.0) return;
    } else if (mode == 2) {
        // Reverse the path.
        xyzout.close();
        xyzin.open(xyz_fn);
        std::vector<std::string> store_hof;
        std::vector<std::vector<std::string>> store_path;
        {
            std::string dummy, buf;
            int i = 0;
            while (i < 100000) {
                if (!std::getline(xyzin, dummy)) break;
                if (!std::getline(xyzin, buf)) break;
                store_hof.push_back(buf);
                std::vector<std::string> frame;
                for (int k = 0; k < nl_atoms; ++k) {
                    if (!std::getline(xyzin, buf)) break;
                    frame.push_back(buf);
                }
                store_path.push_back(frame);
                ++i;
            }
        }
        xyzin.close();
        xyzout.open(xyz_fn, std::ios::out | std::ios::trunc);
        int nframes = (int)store_hof.size();
        npt = 0;
        for (int i = nframes - 1; i >= 1; --i) {  // exclude first (common) point
            xyzout << std::setw(6) << nl_atoms << " \n";
            size_t kpos = store_hof[i].find("HEAT");
            ++npt;
            xyzout << std::setw(8) << store_hof[i].substr(0, 8) << " "
                   << std::setw(4) << npt << " " << store_hof[i].substr(kpos)
                   << "\n";
            for (int k = 0; k < nl_atoms; ++k)
                xyzout << store_path[i][k] << "\n";
        }
        imodel = 0;
    }
}
