// ============================================================================
// 【中文导读】esp.h —— 偶极/极化/静电接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// esp.h — C++ interface for MOPAC 2016 "esp.F90" (ESP charge fitting).
#pragma once

void esp();
void pdgrid();
void surfac();
void potcal();
void elesn();
void espfit();
void naicas();
void ovlp(int ic);
void setup3();
void fsub(int n, double x, double& fval);
void naicap();
void setup_esp(int mode);
void getattrib(double xmin[4], double xmax[4]);
void espplane(int iplane, const double xmin[4], const double step[4],
              int ngridpts2, int ngridpts3);
double dex2(int m);
