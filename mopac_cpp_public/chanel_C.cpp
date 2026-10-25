// ============================================================================
// 【中文导读】chanel_C.cpp —— I/O 通道与文件名表（对应 chanel_C.F90）
// ----------------------------------------------------------------------------
// 定义 Fortran 版 unit 号（iw/ir/iden/ires/iarc/ixyz/...）及各类文件名
// （input_fn/output_fn/restart_fn/density_fn/archive_fn/cosmo_fn/...）与
// ioda 表。C++ 用 std::string/vector 替代 Fortran character 数组，数值语义一致。
// ============================================================================

// ============================================================================
// 【中文导读】chanel_C.cpp —— I/O 通道与文件名表（对应 chanel_C.F90）
// ----------------------------------------------------------------------------
// 定义 Fortran 版 unit 号（iw/ir/iden/ires/iarc/ixyz/...）及各类文件名
// （input_fn/output_fn/restart_fn/density_fn/archive_fn/cosmo_fn/...）与
// ioda 表。C++ 用 std::string/vector 替代 Fortran character 数组，数值语义一致。
// ============================================================================

// chanel_C.cpp — storage.
#include "chanel_C.h"
#include <vector>
namespace chanel_C {
int iw = 6;
int ir = 5;
int imep = 0; int iden = 4;
int igpt = 0; int iw0 = 0; int ires = 9; int iarc=8, ixyz=10;
bool log = false;
int ilog = 7;
std::string input_fn, output_fn, restart_fn, density_fn, log_fn, end_fn;
std::string archive_fn, brillouin_fn, esp_fn, ump_fn, mep_fn, pol_fn;
std::string gpt_fn, esr_fn, xyz_fn, syb_fn, cosmo_fn;
std::string job_fn;
int iscr = 12;
int lbfgs_it = 13;
std::FILE* out_fp = nullptr;   // nullptr => stdout
int iump = 0;
std::vector<int> ioda(146,0);  // F90 chanel_C: dimension(145)
std::vector<int> ifilen(146,0);  // F90 chanel_C: dimension(145)
int irecst = 0;
int idaf = 17;
int irecln = 1023;
}
