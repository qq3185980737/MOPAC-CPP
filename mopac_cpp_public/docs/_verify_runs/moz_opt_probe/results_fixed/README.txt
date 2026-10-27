results_fixed — MOZYME 几何优化修复后的验证记录（2026-10-26，§100）
=================================================================

修复内容
  1) dcart 的 MOZYME 分支密度块索引改为 0-based（与 dhc() 的 ABI 一致）；
  2) hcore_for_MOZYME 恢复 addnucz(phinet, qscnet, qdenet) 调用（COSMO 梯度用）。

结果
  * MOZYME 优化（气体/COSMO）首轮梯度、最终 HOF 与官方一致：
      水        3.818 / 7.768（COSMO），HOF -57.79909 / -65.21180
      苯        17.244 / 20.182（COSMO），HOF 22.95671 / 20.95883
      萘        119.087 / 121.751，HOF 39.71539 / 36.32591
      环戊二烯  210.033 / 232.545，HOF 32.29160 / 68.32543
  * 既有 14 项回归全部 PASS（THERMO 10.007 / FORCE 10.120 不变）；
  * COSMO 套件 12 项中 10 项完全一致（water/meoh MOZYME+COSMO 优化由 FAIL 转 PASS）。
  * 仍剩：乙酸根 1e-5；甲醇普通 COSMO 优化 1e-3（EF 优化器精度）；
    C12H26 大体系 EF 收敛差异（非 MOZYME，非 MOZYME 同样差）。

目录
  mozfix_probe2.log   修复后 old_probe.py 运行记录（水/苯/C12 气体+COSMO）
  gates_mozfix.log    14 项回归
  cosmo_mozfix.log    COSMO 套件
  results_fixed/<case>/{cpp,ref}/  .mop 与 .out（水/苯/萘/环戊二烯）
