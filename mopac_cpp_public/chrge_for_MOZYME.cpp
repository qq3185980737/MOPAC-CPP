// ============================================================================
// 【中文导读】chrge_for_MOZYME.cpp —— 结构/成键/电荷处理
// ----------------------------------------------------------------------------
// 对应 MOPAC 2016 chrge_for_MOZYME.F90 的 C++ 移植（或同源子程序）。
// 约定：保持 Fortran 1-based 索引、打包下三角/列主序布局；具体实现与
// 关键字分支见文件内函数及原有英文注释；缺失/占位处已记录在
// docs/TRANSLATION_GAPS.md。
// ============================================================================

// chrge_for_MOZYME.cpp
#include "chrge_for_MOZYME.h"
extern int ijbo(int i, int j);
void chrge_for_MOZYME(const double* p, double* q, int numat, const int* iorbs) {
    for (int i=1;i<=numat;++i) {
        int ii=ijbo(i,i);
        double sum=0.0;
        for (int j=1;j<=iorbs[i];++j) { ii+=j; sum+=p[ii]; }
        q[i]=sum;
    }
}
