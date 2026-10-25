// ============================================================================
// 【中文导读】polar.h —— 偶极/极化/静电接口（C++ 头文件）
// ----------------------------------------------------------------------------
// 声明对应的模块全局量/函数/常量，与同名 .cpp 配对使用；保持 MOPAC
// 的 1-based 索引与打包下三角/列主序约定。实现细节与移植缺口见
// docs/TRANSLATION_GAPS.md 及对应 .cpp 的中文导读。
// ============================================================================

// polar.h — C++ translation.
#pragma once
#include <vector>
#include <string>
void polar();
void zerom(std::vector<std::vector<double>>& x, int m);
void tf(std::vector<std::vector<double>>& ua, std::vector<std::vector<double>>& ga,
        std::vector<std::vector<double>>& ub, std::vector<std::vector<double>>& gb,
        std::vector<std::vector<double>>& t, int norbs);
void transf(std::vector<std::vector<double>>& f, std::vector<std::vector<double>>& g,
            std::vector<std::vector<double>>& c, int norb);
double trsub(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
             std::vector<std::vector<double>>& ur, int l1, int lm, int ndim);
double trudgu(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
              std::vector<std::vector<double>>& ur, int l1, int lm, int ndim);
double trugdu(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
              std::vector<std::vector<double>>& ur, int l1, int lm, int ndim);
double trugud(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
              std::vector<std::vector<double>>& ur, int l1, int lm, int ndim);
double wrdkey(const std::string& keywrd_, const std::string& key, int nk,
              const std::string& refkey, int nr, double def);
double aval(std::vector<std::vector<double>>& h, std::vector<std::vector<double>>& d, int norbs);
double pol_vol(double average);