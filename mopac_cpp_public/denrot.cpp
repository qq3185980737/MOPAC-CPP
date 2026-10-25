// denrot.cpp — C++ translation of MOPAC 2016 "denrot.F90".
// Prints the density matrix in the local (s-sigma, p-sigma, p-pi) basis.
#include "denrot.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "coe.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "gmetry.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;

extern void denrot_for_MOZYME();

namespace {
// irot(5,35): flat data from denrot.F90 (column-major).
const int kIrotFlat[175] = {
    1,1,1,3,3,  2,2,2,4,3,  3,2,2,2,3,  4,2,2,3,3,  2,3,2,4,2,
    3,3,2,2,2,  4,3,2,3,2,  2,4,2,4,4,  3,4,2,2,4,  4,4,2,3,4,
    5,5,3,1,5,  6,5,3,4,3,  7,5,3,3,3,  8,5,3,2,3,  9,5,3,5,3,
    5,6,3,1,2,  6,6,3,4,2,  7,7,3,3,2,  8,6,3,2,2,  9,6,3,5,2,
    5,7,3,1,4,  6,7,3,4,4,  7,7,3,3,4,  8,7,3,2,4,  9,7,3,5,4,
    5,8,3,1,1,  6,8,3,4,1,  7,8,3,3,1,  8,8,3,2,1,  9,8,3,5,1,
    5,9,3,1,5,  6,9,3,4,5,  7,9,3,3,5,  8,9,3,2,5,  9,9,3,5,5};
const int kIsp[10] = {0, 1, 2, 3, 3, 4, 5, 5, 6, 6};
const char* kAtorbs[9] = {"S-SIGMA", "P-SIGMA", "  P-PI ", "  P-PI ",
                          "D-SIGMA", "  D-PI ", "  D-PI ", " D-DELL",
                          " D-DELL"};
}  // namespace

void denrot() {
    if (mozyme) {
        denrot_for_MOZYME();
        return;
    }
    std::vector<double> b(mpack + 1, 0.0);
    gmetry(geo, coord);
    std::fprintf(stderr, "DENROT after gmetry iprt0\n");

    int iprt = 0;
    std::vector<int> natom(mpack + 1, 0);
    std::vector<std::string> itext(mpack + 1, " "), jtext(mpack + 1, " ");

    for (int i = 1; i <= numat; ++i) {
        int i_f = nfirst[i];
        int i_l = nlast[i];
        int ipq = i_l - i_f - 1;
        int ii = ipq + 2;
        if (ii == 0) continue;
        for (int i1 = 1; i1 <= ii; ++i1) {
            int j1 = iprt + kIsp[i1];
            itext[j1] = kAtorbs[i1 - 1];
            jtext[j1] = elemnt[nat[i]];
            natom[j1] = i;
        }
        iprt += kIsp[ii];  // F90: iprt = j1 = iprt + isp(ii)
        std::fprintf(stderr, "CXXDEN atom %d ii %d ipq %d iprt %d\n", i, ii, ipq, iprt);
        if (ipq != 2) ipq = std::min(std::max(ipq, 1), 3);

        for (int j = 1; j <= i; ++j) {
            std::fprintf(stderr, "DENROT i=%d j=%d\n", i, j);
            double delx = coord[0][j] - coord[0][i];
            double dely = coord[1][j] - coord[1][i];
            double delz = coord[2][j] - coord[2][i];
            if (id != 0) {
                double rmin = 100.0, delxm = 0, delym = 0, delzm = 0;
                for (int mm = -l1u; mm <= l1u; ++mm)
                    for (int nn = -l2u; nn <= l2u; ++nn)
                        for (int kk = -l3u; kk <= l3u; ++kk) {
                            double dx = tvec[0][0] * mm + tvec[0][1] * nn + tvec[0][2] * kk + delx;
                            double dy = tvec[1][0] * mm + tvec[1][1] * nn + tvec[1][2] * kk + dely;
                            double dz = tvec[2][0] * mm + tvec[2][1] * nn + tvec[2][2] * kk + delz;
                            double rr = dx * dx + dy * dy + dz * dz;
                            if (rr >= rmin) continue;
                            rmin = rr; delxm = dx; delym = dy; delzm = dz;
                        }
                delx = delxm; dely = delym; delz = delzm;
            }
            int j_f = nfirst[j];
            int j_l = nlast[j];
            int jpq = j_l - j_f - 1;
            int jj = jpq + 2;
            if (jj == 0) continue;
            if (jpq != 2) jpq = std::min(std::max(jpq, 1), 3);

            std::vector<std::vector<double>> pab(10, std::vector<double>(10, 0.0));  // F90 pab(9,9)
            int kk = 0;
            for (int k = i_f; k <= i_l; ++k) {
                ++kk;
                for (int l = 1; l <= jj; ++l) pab[kk][l] = p[j_f + (k * (k - 1)) / 2 + l - 1];
            }

            double c[75] = {};
            double rdist = 0.0;
            std::fprintf(stderr, "DENROT before coe ii=%d jj=%d ipq=%d jpq=%d\n", ii, jj, ipq, jpq);
            coe(delx, dely, delz, ipq, jpq, c, rdist);
            std::fprintf(stderr, "DENROT after coe\n");
            std::vector<std::vector<double>> arot(10, std::vector<double>(10, 0.0));
            for (int i1 = 1; i1 <= 35; ++i1) {
                int r1 = kIrotFlat[(i1 - 1) * 5 + 0];
                int r2 = kIrotFlat[(i1 - 1) * 5 + 1];
                int c1 = kIrotFlat[(i1 - 1) * 5 + 2];
                int c2 = kIrotFlat[(i1 - 1) * 5 + 3];
                int c3 = kIrotFlat[(i1 - 1) * 5 + 4];
                // F90 c(3,5,5) column-major: c(i)+3*(j-1+5*(k-1)), 1-based i,j,k
                arot[r1][r2] = c[(c1 - 1) + 3 * ((c2 - 1) + 5 * (c3 - 1))];
            }
            int l1 = kIsp[ii];
            int l2 = kIsp[jj];
            std::vector<std::vector<double>> vect(10, std::vector<double>(10, -1.0));
            if (l2 > 0)
                for (int a = 1; a <= l1; ++a)
                    for (int bb = 1; bb <= l2; ++bb) vect[a][bb] = 0.0;
            if (i != j) {
                int ij = std::max(ii, jj);
                for (int i1 = 1; i1 <= ii; ++i1)
                    for (int j1 = 1; j1 <= jj; ++j1) {
                        double sum = 0.0;
                        for (int a = 1; a <= ij; ++a)
                            for (int bb = 1; bb <= ij; ++bb)
                                sum += arot[a][i1] * pab[a][bb] * arot[bb][j1];
                        vect[kIsp[i1]][kIsp[j1]] += sum * sum;
                    }
            }
            kk = 0;
            for (int i1 = i_f; i1 <= i_l; ++i1) {
                ++kk;
                int l = 0;
                for (int j1 = j_f; j1 <= j_l; ++j1) {
                    ++l;
                    if (j1 > i1) continue;
                    b[j1 + (i1 * (i1 - 1)) / 2] = vect[kk][l];
                }
            }
        }
    }

    std::fprintf(stderr, "DENROT after loops iprt=%d mpack=%d\n", iprt, mpack);
    {
        std::fprintf(stderr, "CXXDEN b");
        for (int k = 1; k <= std::min(9, mpack); ++k) std::fprintf(stderr, " %.10f", b[k]);
        std::fprintf(stderr, "\n");
    }
    // Remove dead space (entries <= -0.1).
    int lcount = 0;
    for (int i = 1; i <= mpack; ++i) {
        if (b[i] <= -0.1) continue;
        ++lcount;
        b[lcount] = b[i];
    }

    // Put atomic orbital valencies onto the diagonal.
    for (int i = 1; i <= iprt; ++i) {
        double sum = 0.0;
        int ii = (i * (i - 1)) / 2;
        for (int j = 1; j <= i; ++j) sum += b[j + ii];
        for (int j = i + 1; j <= iprt; ++j) sum += b[(j * (j - 1)) / 2 + i];
        b[(i * (i + 1)) / 2] = sum;
    }

    // Print the rotated density matrix.
    int limit = (iprt * (iprt + 1)) / 2;
    int kk = 8, na = 1;
    while (true) {
        int ll = 0;
        int m = std::min(iprt + 1 - na, 6);
        int ma = 2 * m + 1;
        m = na + m - 1;
        std::printf("\n%16s", "");
        for (int i = na; i <= m; ++i) std::printf(" %7s   ", itext[i].c_str());
        std::printf("\n%15s", "");
        for (int i = na; i <= m; ++i) std::printf("  %2s%3d    ", jtext[i].c_str(), natom[i]);
        std::printf("\n");
        for (int k = 1; k <= ma; ++k) std::printf("%6s", "------");
        std::printf("\n");
        int l_hi = 0;
        for (int i = na; i <= iprt; ++i) {
            ++ll;
            int k0 = (i * (i - 1)) / 2;
            l_hi = std::min(k0 + m, k0 + i);
            int k_start = k0 + na;
            if (kk + ll > 50) {
                std::printf("\n%17s", "");
                for (int n = na; n <= m; ++n) std::printf(" %7s   ", itext[n].c_str());
                std::printf("\n%17s", "");
                for (int n = na; n <= m; ++n) std::printf("  %2s%3d    ", jtext[n].c_str(), natom[n]);
                std::printf("\n");
                for (int n = 1; n <= ma; ++n) std::printf("%6s", "------");
                std::printf("\n");
                kk = 4; ll = 0;
            }
            std::printf(" %7s %2s%3d", itext[i].c_str(), jtext[i].c_str(), natom[i]);
            for (int n = k_start; n <= l_hi; ++n) std::printf("%11.6f", b[n]);
            std::printf("\n");
        }
        if (l_hi >= limit) break;
        kk = kk + ll + 4;
        na = m + 1;
        if (kk + iprt + 1 - na > 50) kk = 4;
    }
}
