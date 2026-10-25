moz_opt_probe — MOZYME 几何优化可用性实测（2026-10-26）
=======================================================

结论：**不可用**。MOZYME 几何优化（气体与 COSMO）在 C++ 端从第 1 轮起梯度就
错误约 2 个数量级，优化器无法推进；能量（首轮 HOF）与官方完全一致，说明问题在
梯度层（dcart 的 MOZYME 分支），不在优化器、也不在本轮 COSMO 修复。

运行方法
--------
    python docs\_verify_runs\moz_opt_probe\run_moz_opt_probe.py --all

（需要 bin\mopac_cpp.exe 与 bin\mopac_ref.exe，见仓库说明；每个用例两种程序
各跑一次，输出写入 results\<case>\{cpp,ref}\。）

用例与结果（首轮梯度 g1 / 最终 HOF / 最终梯度范数）
--------------------------------------------------
  水 默认           C++ g1=721.4  HOF −57.79108 GN 722.6 | ref g1=3.82  HOF −57.79909 GN 0.78
  水 BFGS           C++ 0 轮      HOF −57.78992          | ref 3 轮   HOF −57.79922
  水 L-BFGS         C++ 同默认（关键字未生效）           | ref 报错（关键字不支持）
  水 EF             同默认                               | 同默认
  水 COSMO          C++ g1=684.1  HOF −65.17074 GN 693.8 | ref g1=7.77  HOF −65.21180 GN 0.46
  苯 默认           C++ g1=6284.3 HOF 23.02282  GN 6284 | ref g1=17.24 HOF 22.95671 GN 0.58
  苯 BFGS           C++ 0 轮      HOF 23.16493          | ref HOF 22.96272
  苯 COSMO          C++ g1=6220.8 HOF 21.05900 GN 6220 | ref g1=20.18 HOF 20.95883 GN 0.67
  C12H26 默认       C++ g1=5513.4 HOF 1252.10  GN 5594 | ref g1=510.9 HOF 162.75 GN 0.99
  C12H26 BFGS       C++ 0 轮 FAILED  HOF 1231.60        | ref 63 轮 HOF 171.14
  C12H26 COSMO      C++ g1=5470.6 HOF 1234.42           | ref g1=514.1 HOF 98.40

关键点
------
* 首轮能量（HOF）C++ 与官方完全一致（水 −57.78972），SCF/能量路径正常；
* 首轮梯度即错 2 个数量级（水 721 vs 3.8；苯 6284 vs 17.2；C12 5513 vs 511），
  与分子、优化器、是否 COSMO 无关 → 梯度层问题；
* C++ 优化 3 轮后以巨大梯度"结束"，HOF 停在起点附近；BFGS 甚至 0 轮；
* `L-BFGS` 关键字官方不支持（报错），C++ 静默按默认处理；
* `DCART` 打印的笛卡尔导数同样偏大 2 个数量级（水：O 266/−385/−252 vs
  0.05/0.0003/−3.10），定位到 dcart 的 MOZYME 分支。

下一步
------
为 dcart 重建 Fortran dump 对照（复用 docs/_verify_runs 中 ifx 参考构建方法），
逐项比对 MOZYME 分支的 padi/pdi/pbdi、dhc 的 aa/ee、delsta 的 dstat 后修复。

目录内容
--------
  run_moz_opt_probe.py   测试脚本（15 个用例 × cpp/ref）
  results/<case>/{cpp,ref}/  .mop 输入与 .out/stderr 记录
