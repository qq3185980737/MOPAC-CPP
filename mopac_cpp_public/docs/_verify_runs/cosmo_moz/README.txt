cosmo_moz — MOZYME / 普通 COSMO 回归记录（2026-10-26）
======================================================

用途
----
验证 linear_cosmo / mkl_bits / cosmo(diegrd) 修复后 COSMO 相关数值：
  * MOZYME + COSMO 1SCF；
  * 普通（非 MOZYME）COSMO 1SCF 与几何优化；
并记录仍存在的已知问题。

运行方法
--------
1) 准备可执行文件（bin/，*.exe 已被 .gitignore 忽略）：
   x64\Release\mopac.exe                      -> bin\mopac_cpp.exe
   D:\test_rep\mopac_fortran_ref\mopac_ref.exe -> bin\mopac_ref.exe
2) 在仓库根目录运行：
   python docs\_verify_runs\cosmo_moz\run_cosmo_moz.py --all

结果（2026-10-26，修复后）
--------------------------
MOZYME+COSMO 1SCF（PM7 MOZYME 1SCF EPS=78.4）：
  水            HOF −65.16392  DIEL −0.36147   与官方一致
  甲醇          HOF −53.32319  DIEL −0.29319   与官方一致
  苯            HOF  21.05947  DIEL −0.09229   与官方一致
  环戊二烯      HOF 105.83580  DIEL −0.77391   与官方一致
  萘            HOF  39.79622  DIEL −0.15878   与官方一致
  乙酸根(CHARGE=-1) HOF 47.18759 / −3.37822    官方 47.18758（差 1e-5，待查）

非 MOZYME COSMO 1SCF（PM7 1SCF EPS=78.4）：
  水            HOF −65.16522  DIEL −0.35863   与官方逐位一致
  甲醇          HOF −53.32397  DIEL −0.29253   与官方逐位一致

非 MOZYME COSMO 几何优化（PM7 EPS=78.4，EF 优化器）：
  水            首轮梯度 8.388，HOF −65.21143 / DIEL −0.36538   与官方完全一致
  甲醇          首轮梯度 32.005 与官方一致，收敛 HOF −54.01839 vs −54.01739
                （差 1e-3，EF 优化器已知精度，docs §72-82）

已知问题
--------
* MOZYME 几何优化（含 COSMO）本身不可用：首轮梯度 684/3451 vs 官方 7.77/31.88，
  属既有 MOZYME 梯度链路问题（docs §71 已记录；气体 MOZYME 优化同样不收敛），
  与本轮 diegrd/COSMO 修复无关。COSMO 优化崩溃（0xC0000005）已在本轮修复。
* 乙酸根 HOF 1e-5 残差（DIELECTRIC 已一致）。

目录内容
--------
  run_cosmo_moz.py   回归脚本（6 个 1SCF + 2 个优化用例）
  decks/*.mop        输入 deck
  results/           本次 cpp / ref 输出记录
  bin/               放置 mopac_cpp.exe / mopac_ref.exe
