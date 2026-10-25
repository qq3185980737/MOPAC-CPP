// ============================================================================
// 【中文导读】lbfgs_line_search.h —— MOPAC 2016 模块接口声明（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

#include <limits>
void dcstep(double& stx,double& fx,double& dx,double& sty,double& fy,double& dy,double& stp,double fp,double dp,bool& brackt,double stpmin,double stpmax);
void dcsrch(double f,double g,double& stp,double ftol,double gtol,double xtol,double stpmin,double stpmax,std::string& task,int* isave,double* dsave);
double dpmeps();