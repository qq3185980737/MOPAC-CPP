// ============================================================================
// 【中文导读】derivs_C.cpp —— 导数模块的全局存储（derivs_C.F90）
// ----------------------------------------------------------------------------
// 存放解析导数（deri1/deri2/deri21-23/dernvo）共用的超向量：b/ab/fb_ci/
// fmat/hmat/wmat 等；全部 1-based 打包/列主序，尺寸按 meci/molkst 参数确定。
// ============================================================================

// derivs_C.cpp — storage for derivs_C.
#include "derivs_C.h"
namespace derivs_C {
std::vector<double> b;
std::vector<double> ab;
std::vector<double> fb_ci;
std::vector<double> fmat;
std::vector<double> hmat;
std::vector<double> wmat;
std::vector<double> aidref;
std::vector<double> work2;
}
