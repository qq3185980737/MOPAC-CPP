# MOZYME+MD 能量守恒修复与探究记录

> 研究主体：168 原子氢氧热点点火体系（42 O2 + 42 H2，4 热点 O2，BOX=175，无 VELO）
> 主线目标：自研 C++ MOPAC（MOZYME + MD）总能量守恒。量化门槛未最终拍板；
> 实测口径：`Total Energy (Potential + Kinetic) =` 序列首末差 ×627.509 ÷（帧数×dt fs），单位 kcal/mol/fs。

## 1. 现象与根因（三环泄漏链）

首跑（1000 步/300 fs，默认 SCFCRT=1e-2）总能量漂移 **-3.84 kcal/mol/fs**，成水 607/1000 帧。
SCFCRT=1e-5 试跑漂移 -6.69（更差）——SCFCRT 关键词在 MOZYME 热启动下"无效"。

探针（MOPAC_MOZ_SCFDBG，3 步 2013 条，全部 step=1）锁定根因：每帧 SCF 只 3–4 次迭代即 okscf=1
（ovmax≈0.09–0.13、ediff≈0.6–1.5 eV≈1 kcal 量级），走的是 isitsc `iemin/iemax`"能量平台"判据
（相对容差 `0.1*(emin-escf)`），完全绕开 selcon 判据。iscf=1 只是收敛类别标签，非迭代次数。

根因是三条互相独立的跨帧状态泄漏：

1. **md_driver 不递增 step_num**（递增点原本仅在 run_mopac.cpp:245、set_up_RAPID.cpp:48、
   big_swap.cpp:523/532）→ iter_for_MOZYME.cpp:234 `if (icalcn != step_num)` 的每步重置块
   （eold/ovmax/scfcri/idiagg/tiny）只在 MD 首帧执行一次，之后每帧热启动。
2. **compfg.cpp:166 `emin=0.0` 只在首帧（l1==0）执行；:327 每帧把 emin 更新为历史最低**
   → emin 跨帧保留后，isitsc 平台判据容差 `0.1*(emin-escf)≈34 kcal` 巨大，热启动快速"平台收敛"。
3. **isitsc.cpp 内部 static scf1/escf0[6] 跨帧保留**，能量平台历史被上一帧"借用"。

## 2. 修复落点（5 处，已重建验证）

| 位置 | 改动 | 作用 |
|---|---|---|
| isitsc.h | 声明 `isitsc_reset()` | 接口 |
| isitsc.cpp | static scf1/escf0[6] 提升文件级 + 实现 `isitsc_reset()` | 每帧清零平台历史 |
| iter_for_MOZYME.cpp 重置块 | 追加 `iemin=0; iemax=0; emin=0.0; isitsc_reset();` | emin=0 使 isitsc 首调用 okscf=false（`if(emin==0.0){okscf=false;return;}`），平台判据整体禁用，SCF 只能靠密度残差判据或撞 itrmax |
| md_driver.cpp 主循环首行 | `step_num = step_num + 1;` | 每帧触发重置块 + 密度重建，LMO 系数热启动保留 |
| 探针块 | MOPAC_MOZ_SCFDBG fprintf 分支移除 | 交付态干净 |

判据语义（isitsc.cpp）：`okscf = (ovmax<selcon*5 && |ediff|<selcon && scf1) || (niter>itrmax)`。

## 3. 实测守恒序列（全部本机实测）

| 配置 | 帧/物理时间 | 漂移 kcal/fs | 每帧迭代 | wall | 备注 |
|---|---|---|---|---|---|
| 首跑（泄漏，1e-2） | 1000/300 fs | **-3.84** | ~4（平台） | 128 s | 成水 607/1000、产率 2.3%、T 480→5491 K |
| SCFCRT=1e-5（证伪） | 300/90 fs | -6.69 | 2001 撞 itrmax | 52 s | ovmax 停 1.37e-3（MOZYME 截断下限），该方向封死 |
| 修复 + 1e-2 | 300/90 fs | -1.30 | 31–39 | 107 s | E0=5.388200/Elast=5.202300 |
| 修复 + 1e-3 | 100/30 fs | -0.57 | ~45 | 135 s | |
| **修复 + 5e-4** | 300/90 fs | **-0.33** | ~65 | 586 s | E0=5.386400/Elast=5.339200；T 末帧 4835 K、峰值 5598 K；成水 437/300 帧、末段 1.64 水/帧、产率 3.9% |
| **修复 + 5e-4 + DT=0.1** | 900/300 fs | **-0.144** | ~65 | 722 s | E0=5.3895/Elast=5.3688；T 末帧 4689 K；总漂移 -12.99 kcal=总能量 0.4% |

结论：漂移 ≈ ∝ sqrt(selcon)；缩小步长对 SCF 噪声注入 ∝ dt（0.3→0.1 fs 改善 2.3×）。

## 4. 方法探究结论（用户询问"增加迭代次数/缩小步长"）

- **增加迭代次数**：迭代数由 SCFCRT 决定，非独立参数；selcon=1e-5 时 2000 次迭代 ovmax 仍
  停 1.37e-3（MOZYME 局部化截断的密度残差物理下限），放开 itrmax 无增益——已实测封死。
- **缩小步长**：md_driver 已是 dt^4 标度高阶传播器（3 力评估/步，md_driver.cpp:479），积分误差
  占比小；但 SCF 噪声注入 ∝ dt（同物理时间动能漂移 ∝ δF^2·T·dt），0.3→0.1 fs 实测改善 2.3×。
- **判据改造**（ovmax 放宽固定 5e-3 + ediff 收紧 1e-4）：偏离官方 isitsc 语义（selcon 是密度∧能量
  同尺度统一标尺）；且解析梯度力噪声 ∝ 密度残差 δP，ovmax 放宽可能反而增大力噪声，净效果不可
  先验断定——**未获授权，已搁置**。
- **接受水平**：DT=0.1 + SCFCRT=5e-4 → -0.144 kcal/fs（总能量 0.4%），作为当前主线交付水平。

## 5. 崩溃排查（dt=0.1 首跑 632/900 帧硬崩）

- 现象：exit -1、无错误消息、无 SCF restart/jump 痕迹、out 停在坐标块中途、stderr 尾部正常 deriv 行。
- 判定：50 步 smoke 正常（exit=0、0 restart）；重跑完整 900 帧 exit=0、内存恒定 45 MB（排除泄漏
  与 dt 相关确定性 bug）→ **偶发崩溃**，未再复现。
- 隐患记录：md_driver.cpp:795 SCF-restart 的 `dt *= 0.5` 无下界（:905 jump 重试有 0.01 fs 下界），
  未触发，但属潜在风险点。

## 6. 待办 / 影响面

- **非 MD 位级回归（2026-11-03 已完成）**：发布语料 `regression/` 349 输入不含 MOZYME 关键词；
  用交付态 exe（commit a3f79e7）重跑全部 327 个有 Fortran .arc 对照的输入（OMP/MKL 单线程）：
  **326/327 与 out_ftn 位级一致**，唯一差异 `_val_ef`（EF 优化端点 ULP：-4.26000 vs -4.25995），
  0 崩溃/0 缺失 —— 修复对非 MOZYME 路径零破坏。基线（记录输出）为 320/327 一致、7 个 EF ULP；
  重跑后 6 个 `_mx_opt_*` EF 优化转一致，反映发布快照之后 EF 链修复（fock2 hl1 结合序对齐等）
  已并入当前 exe。MOZYME 专项：6 个 MOZYME+COSMO 1SCF 与官方逐位一致（cosmo_moz 记录不变）。
- **MOZYME 几何优化意外反转（2026-11-03 实测）**：moz_opt_probe 曾记录"C++ 端 MOZYME 优化不可用"
  （首轮梯度错 2 数量级）。修复后复测：SCF 平台假收敛解除 → 首轮梯度与官方一致（水 3.818 vs
  3.82、苯 17.244 vs 17.24、C12 510.9 vs 510.9），水/苯优化收敛 HOF 与官方逐位一致
  （-57.79909 / 22.95671），C12 177 轮收敛 GN 0.968、HOF 163.55（官方 162.75，EF 轨迹敏感性级）。
  **即三环泄漏修复同时解除了 MOZYME 梯度错误的根因**；moz_opt_probe/README.txt 已追加复测记录。
- 1000 步长程守恒未跑（当前 900 步/300 fs 已闭环）。
- 修复代码 4 文件已提交；产物/脚本在 D:\mopac_md_work\mozmd\（moz168_fix5.out、moz168_fix_dt.out 等）。
