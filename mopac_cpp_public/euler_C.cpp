// ============================================================================
// 【中文导读】euler_C.cpp —— 周期胞整数下界（euler_C.F90）
// ----------------------------------------------------------------------------
// 定义 l1l/l2l/l3l（周期胞在三个方向的整数下界），与 euler_C.h 的其余成员
// 一起供周期性计算（最小镜像/晶胞求和）使用。
// ============================================================================

// euler_C.cpp — C++ translation of MOPAC 2016 "euler_C.F90".
// Defines l1l/l2l/l3l (periodic-cell lower bounds); see euler_C.h for the
// split of the remaining module members.
#include "euler_C.h"

namespace euler_C {
int l1l = 0, l2l = 0, l3l = 0;
}  // namespace euler_C
