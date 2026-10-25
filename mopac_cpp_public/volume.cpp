// ============================================================================
// 【中文导读】volume.cpp —— 网格/体积/静电辅助
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 volume.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// volume.cpp
#include "volume.h"
#include <cmath>
double volume(double* v, int ndim) {
    double a=std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    if (ndim==1) return a;
    double b=std::sqrt(v[3]*v[3]+v[4]*v[4]+v[5]*v[5]);
    double dx=v[0]-v[3],dy=v[1]-v[4],dz=v[2]-v[5];
    double g=std::sqrt(dx*dx+dy*dy+dz*dz);
    double c=(a*a+b*b-g*g)/(2.0*a*b);
    double sing=std::sqrt(1.0-c*c);
    if (ndim==2) return a*b*sing;
    return std::abs((v[1]*v[5]-v[2]*v[4])*v[6]+(v[2]*v[3]-v[0]*v[5])*v[7]+(v[0]*v[4]-v[1]*v[3])*v[8]);
}
