// ============================================================================
// 【中文导读】elemts_C.cpp —— 元素符号/名称表（elemts_C.F90）
// ----------------------------------------------------------------------------
// 1-based 元素表：elemnt[1..107] 为 MOPAC 固定两字符符号（如 "H "、"CL"），
// cap_elemnt 为大小写混合形式，另含原子量 atmass。表的下标即原子序数，
// 与 nat[] 的取值一一对应；打印/解析元素时必须保持两字符（含尾随空格）。
// ============================================================================

// elemts_C.cpp — storage.
#include "elemts_C.h"
namespace elemts_C {
std::vector<std::string> elemnt = {
    "",  // index 0 placeholder (table is 1-based: elemnt[1]="H ")
    "H ", "HE", "LI", "BE", "B ", "C ", "N ", "O ", "F ", "NE",
    "NA", "MG", "AL", "SI", "P ", "S ", "CL", "AR", "K ", "CA", "SC", "TI",
    "V ", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR",
    "KR", "RB", "SR", "Y ", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG",
    "CD", "IN", "SN", "SB", "TE", "I ", "XE", "CS", "BA", "LA", "CE", "PR",
    "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF",
    "TA", " W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO",
    "AT", "RN", "FR", "RA", "AC", "TH", "PA", " U", "NP", "PU", "AM", "CM",
    "BK", "CF", "XX", "+3", "-3", "CB", "++", " +", "--", " -", "TV"};
std::vector<std::string> cap_elemnt = {
    "",  // index 0 placeholder (1-based)
    "H ", "HE", "LI", "BE", "B ", "C ", "N ", "O ", "F ", "NE",
    "NA", "MG", "AL", "SI", "P ", "S ", "CL", "AR", "K ", "CA", "SC", "TI",
    "V ", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR",
    "KR", "RB", "SR", "Y ", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG",
    "CD", "IN", "SN", "SB", "TE", "I ", "XE", "CS", "BA", "LA", "CE", "PR",
    "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF",
    "TA", " W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO",
    "AT", "RN", "FR", "RA", "AC", "TH", "PA", " U", "NP", "PU", "AM", "CM",
    "BK", "CF", "XX", "+3", "-3", "CB", "++", " +", "--", " -", "TV"};
std::vector<std::string> atom_names(108, "            ");
}  // namespace elemts_C

// Legacy function-style element-symbol accessor used by prtlmo/prtgra/superd/geout.
// z is a 1-based atomic number; the table is stored 0-based.
const char* elemnt(int z) {
    if (z < 1 || z >= (int)elemts_C::elemnt.size()) return "  ";
    return elemts_C::elemnt[z].c_str();
}
