// ============================================================================
// 【中文导读】sort.cpp —— MOPAC 2016 计算/工具模块
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 sort.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// sort.cpp
#include "sort.h"
void sort(float* val, std::complex<float>* vec, int n) {
    for (int i=1;i<=n;++i) {
        float x=1e9f; int k=i;
        for (int j=i;j<=n;++j) if (val[j]<x) { k=j; x=val[j]; }
        for (int j=1;j<=n;++j) {
            std::complex<float> s=vec[(k-1)*n+(j-1)];
            vec[(k-1)*n+(j-1)]=vec[(i-1)*n+(j-1)];
            vec[(i-1)*n+(j-1)]=s;
        }
        val[k]=val[i];
        val[i]=x;
    }
}
