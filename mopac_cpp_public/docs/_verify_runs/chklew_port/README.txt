chklew_port — CHKLEW 移植回归记录（2026-10-25）
================================================

用途
----
验证 chklew.F90 移植后 C++ 端口的数值与官方一致，覆盖：
  * 基线 6 项：水/苯 π 1SCF、THERMO、FORCE、水/苯 MOZYME；
  * π 体系 8 项（chklew F90 289-618 控制流重译后新增）：
    ethene、acetylene、butadiene、pyridine、naphthalene、cyclopentadiene、
    indole、anthracene。

运行方法（Windows / PowerShell，仓库根目录）
--------------------------------------------
1) 准备两个可执行文件（bin/ 目录，*.exe 已被 .gitignore 忽略，需本地放置）：
   x64\Release\mopac.exe                     -> bin\mopac_cpp.exe
   D:\test_rep\mopac_fortran_ref\mopac_ref.exe -> bin\mopac_ref.exe
   （mopac_ref.exe 为官方 MOPAC2016 源码 ifort+MKL 编译的参考程序）
2) 运行：
   python docs\_verify_runs\chklew_port\run_gate.py --all
   脚本对每个 deck 分别用 cpp / ref 各跑一次（OMP_NUM_THREADS=1），
   输出写入 run\<tag>\{cpp,ref}\，并与 baseline.json 比对。
3) 单独的 π 体系官方对照（不依赖 baseline）：
   python docs\_verify_runs\chklew_port\compare_pi_cases.py

当前基线（见 baseline.json）
----------------------------
  水 BONDS 1SCF        −57.79103    与官方一致
  苯 π 体系 1SCF       246.84781    与官方一致
  THERMO ZPE           10.007       官方 10.003（已知 dcart 低精度路径）
  FORCE ZPE            10.120       官方 10.119（同上）
  水 MOZYME 1SCF       −57.78972    与官方一致
  苯 MOZYME 1SCF       23.02726     与官方一致（2026-10-25 控制流重译后修复）
  ethene               14.29865     与官方一致
  acetylene            57.81121     与官方一致
  butadiene            261.60662    与官方一致
  pyridine             37.15222     与官方一致
  naphthalene          43.04634     官方 43.04633（差 1e-5，下游 MOZYME 精度）
  cyclopentadiene      120.81544    与官方一致
  indole               70.32731     与官方一致
  anthracene           67.10911     与官方一致

目录内容
--------
  run_gate.py           完整回归脚本（6 基线 + 8 π 体系）
  compare_pi_cases.py   仅 π 体系，直接与官方 exe 比对
  baseline.json         基线数值（含各用例官方对照说明）
  decks/*.mop           6 个基线 deck + 2 个 GEOCHK 调试分支 deck
  pi_cases/*.mop        8 个 π 体系 deck
  results/              本次 cpp / ref 输出记录，便于逐行比对
  bin/                  放置 mopac_cpp.exe / mopac_ref.exe

GEOCHK 用例
-----------
decks/benzene_mozyme_geochk*.mop 与 pi_cases 中的 deck 加 GEOCHK 关键字，
即可验证 chklew 的 DEBUG/LARGE 数组转储；苯、萘的 chklew 调试区域与官方
逐字节一致（含前导空列与 (20I4) 换行）。

已知遗留
--------
1) 苯/萘以外，MOZYME + 硫体系仍严重偏差（thiophene、methanethiol），但同几何
   的非 MOZYME PM7 路径逐位一致，且 thiophene 的 chklew 调试转储与官方逐字节
   一致 —— 问题在 MOZYME 的 S 处理，不在 chklew。详见 TRANSLATION_GAPS.md §95.4。
2) naphthalene 的 1e-5 kcal/mol 残差在 MOZYME SCF/LMO 数值精度范围。
