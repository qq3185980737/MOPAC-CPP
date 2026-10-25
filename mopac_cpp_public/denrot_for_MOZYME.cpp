// denrot_for_MOZYME.cpp — C++ translation of "denrot_for_MOZYME.F90".
#include "denrot_for_MOZYME.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "coe.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "ijbo.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;
using namespace MOZYME_C;

namespace {
const int kIrot[5][10] = {
    {1, 2, 3, 4, 2, 3, 4, 2, 3, 4},
    {1, 2, 2, 2, 3, 3, 3, 4, 4, 4},
    {1, 2, 2, 2, 3, 3, 3, 4, 4, 4},
    {3, 4, 2, 3, 4, 2, 3, 2, 3, 4},
    {3, 3, 2, 2, 2, 2, 2, 4, 4, 3}};
const int kIsp[5] = {0, 1, 2, 3, 3};
}  // namespace

void denrot_for_MOZYME() {
    if (maxtxt != 0) {
        std::printf("%-*s%-*s%s\n", maxtxt / 2 + 3, "", maxtxt + 4, "ATOM",
                    "s-s Sig  s-p Sig  p-p Sig  p-p Pi");
    } else {
        std::printf("   ATOM    ATOM    S-S (Sig)   S-P (Sig)   P-P (Sig)    P-P (PI)\n");
    }
    std::printf("\n");

    int imax = maxtxt + 10;
    std::vector<std::vector<double>> pab(5, std::vector<double>(5, 0.0));

    for (int i = 1; i <= numat; ++i) {
        int io = std::min(iorbs[i], 4);
        bool first = true;
        int ipq = std::max(1, io - 2);
        for (int j = 1; j <= numat; ++j) {
            if (i == j || ijbo(i, j) < 0) continue;
            int ij = ijbo(i, j);
            int jo = std::min(iorbs[j], 9);
            int jpq = std::max(1, jo - 2);

            double delx, dely, delz;
            if (id == 0) {
                delx = coord[0][j] - coord[0][i];
                dely = coord[1][j] - coord[1][i];
                delz = coord[2][j] - coord[2][i];
            } else {
                double rmin = 1.0e6;
                delx = dely = delz = 0.0;
                for (int is_loc = -l1u; is_loc <= l1u; ++is_loc)
                    for (int js = -l2u; js <= l2u; ++js)
                        for (int ks = -l3u; ks <= l3u; ++ks) {
                            double xj[3];
                            for (int l = 1; l <= 3; ++l)
                                xj[l - 1] = coord[l - 1][j] + common_arrays_C::tvec[l - 1][0] * is_loc +
                                            common_arrays_C::tvec[l - 1][1] * js + common_arrays_C::tvec[l - 1][2] * ks -
                                            coord[l - 1][i];
                            double r = xj[0] * xj[0] + xj[1] * xj[1] + xj[2] * xj[2];
                            if (r < rmin) {
                                rmin = r;
                                delx = xj[0]; dely = xj[1]; delz = xj[2];
                            }
                        }
            }

            double sum = 0.0;
            for (int ii = 1; ii <= io; ++ii)
                for (int jj = 1; jj <= jo; ++jj) {
                    ++ij;
                    sum += p[ij] * p[ij];
                    pab[ii][jj] = p[ij];
                }

            if (sum < 0.1) continue;

            double c[75] = {};
            double rdist = 0.0;
            coe(delx, dely, delz, ipq, jpq, c, rdist);

            double arot[5][5] = {};
            for (int i1 = 1; i1 <= 10; ++i1)
                arot[kIrot[0][i1 - 1]][kIrot[1][i1 - 1]] =
                    c[(kIrot[2][i1 - 1] - 1) + 3 * ((kIrot[3][i1 - 1] - 1) +
                                                    5 * (kIrot[4][i1 - 1] - 1))];

            double vect[5][5] = {};
            if (i != j) {
                int ijmax = std::max(io, jo);
                for (int i1 = 1; i1 <= io; ++i1)
                    for (int j1 = 1; j1 <= jo; ++j1) {
                        double s = 0.0;
                        for (int l1 = 1; l1 <= ijmax; ++l1)
                            for (int l2 = 1; l2 <= ijmax; ++l2)
                                s += arot[l1][i1] * pab[l1][l2] * arot[l2][j1];
                        vect[kIsp[i1]][kIsp[j1]] += s * s;
                    }
            }

            double sigpi[5] = {0, 0, 0, 0, 0};
            int lsp = 1;
            sigpi[1] = vect[1][1];
            if (jo > 1) {
                lsp = 2;
                sigpi[2] = vect[1][2];
                if (io > 1) {
                    lsp = 4;
                    sigpi[2] += vect[2][1];
                    sigpi[3] = vect[2][2];
                    sigpi[4] = vect[3][3];
                }
            } else if (io > 1) {
                lsp = 2;
                sigpi[2] = vect[2][1];
            }

            if (first) {
                if (maxtxt > 0)
                    std::printf("%6d%s(%s) %6d%s(%s)", i, elemnt[nat[i]].c_str(),
                                txtatm[i].substr(0, maxtxt).c_str(), j,
                                elemnt[nat[j]].c_str(), txtatm[j].substr(0, maxtxt).c_str());
                else
                    std::printf("%6d%s %6d%s", i, elemnt[nat[i]].c_str(), j,
                                elemnt[nat[j]].c_str());
                for (int q = 1; q <= lsp; ++q) std::printf("%9.5f", sigpi[q]);
                std::printf("\n");
                first = false;
            } else {
                if (maxtxt > 0)
                    std::printf("%-*s%6d%s(%s)", imax, "", j,
                                elemnt[nat[j]].c_str(),
                                txtatm[j].substr(0, maxtxt).c_str());
                else
                    std::printf("%8s%6d%s", "", j, elemnt[nat[j]].c_str());
                for (int q = 1; q <= lsp; ++q) std::printf("%9.5f", sigpi[q]);
                std::printf("\n");
            }
        }
    }
}
