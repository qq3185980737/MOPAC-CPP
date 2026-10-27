# 《FORMULAS_AND_VERIFICATION.docx》分支验证方法覆盖审计

审计日期：2026-10-26
审计对象：`docs/FORMULAS_AND_VERIFICATION.docx`（98 KB，573 段、37 表，文件时间 2026-10-25 18:33）
对照基线：`mopac_cpp_public/run_mopac.cpp` 的计算分派、后处理分支与模式开关；
`docs/_verify_runs/` 的实测资产；`docs/TRANSLATION_GAPS.md` 的修复记录。

## 0 判定标准

- **A（有具体验证方法）**：文档给出算例/deck、对照量、参考值、容差或判定标准
  （或"修复前/修复后/官方"三列数值记录）。
- **B（仅有状态/说明）**：只写"待验证/部分/未完成/可运行"，没有参考值或可执行步骤。
- **C（未涉及）**：文档中没有对应分支的验证方法（部分连公式章节也没有）。

注意：文档 §E.5 自定了"全部验证"的定义（每一分支/模式 × 方法 × 元素，双程序同 deck
对照输出，差异不超过打印/数值容差），并在 §D.4、§E.4、§H.3、§H.4、§I.4 明确声明
**"不能声称全部验证"**。本审计按该定义逐分支核对"是否已有可执行的验证方法"。

## 1 一句话结论

**不包含所有分支的验证方法。** 文档对 13 个主分派分支都有公式与章节，但只有约
6 个分支具备可执行的验证方法（其中 EF/FLEPO/L-BFGS/FORCE 仅验证了能量或部分输出）；
其余主分支（TS/SADDLE/NLLSQ/SIGMA/LOCATE-TS/REFINE-TS、PATHS/PATHK、GRID、
DRC 轨迹、周期性、MOLDYN）只有"待验证/扩展"状态；另有约 10 类模式/选项
（PKA、RAPID、蛋白/PDB 流程、HTML、FIELD、EXTERNAL、Sparkles、PMEP 细节、
分析输出细节、MOZYME 优化、COSMO 优化、MOZYME+COSMO 等）在文档中完全没有
验证方法。

## 2 覆盖矩阵

### 2.1 主分派分支（run_mopac.cpp 的 13 个分支；文档 §1 全部列出）

| # | 分支 | 文档公式/章节 | 文档中的验证方法 | 判定 | 备注 |
|---|---|---|---|---|---|
| 1 | SADDLE（react1） | §6 | E.1「待验证」 | B | 无 HCN/HNC 等算例或参考值 |
| 2 | STEP1/GRID | §7 | E.1「待验证」 | B | 无网格扫描 deck/参考值 |
| 3 | latom≠0 / PATHS、PATHK | §7 | E.1「待验证」 | B | 无路径扫描 deck/参考值 |
| 4 | FORCE/THERM/DFORCE | §8 | E.1「未完成」；I.7（ZPE 10.118、三频率、修复前后+官方）、I.11（熵 45.0697） | A-部分 | 仅水 PM7；E.1/H.2 未回写，仍写"不可用" |
| 5 | DRC/IRC | §9 | E.3 water_drc「可运行，轨迹待逐点比」；I.10（IRC 路径点/初速） | B/部分 | 无逐点轨迹参考值；IRC 1e-3~1e-2 差未定容差 |
| 6 | MOLDYN | §9 | E.1「待验证」；H.2「扩展（官方不识别）」 | B | 设计上无法与官方对照，只有"扩展"结论 |
| 7 | NLLSQ / SIGMA | §6 | E.1「待验证」 | B | 无算例/参考值 |
| 8 | 1SCF / 0SCF | §2 | D.2（37 组方法×分子全同）、E.3、G.2/G.3、A/B/F（H₂ 链） | A | 文档中最完整的一类 |
| 9 | LOCATE-TS / REFINE-TS | §6 | E.1「待验证」 | B | 无蛋白/大体系算例 |
| 10 | DFP / FLEPO / BFGS | §4 | E.1「部分」；I.6（水 FLEPO −57.79968、DFP −57.79966） | A-部分 | 仅能量；无几何对照/收敛轨迹 |
| 11 | TS（EF P-RFO）/虚频 | §6 | E.1「待验证」 | B | 仅 §6 描述"用 FORCE 验证恰有一个虚频"，无算例 |
| 12 | L-BFGS-B（GEO_REF） | §5 | E.1「待验证」；I.9（−57.79962 vs −57.79957） | A-部分 | 仅水能量；H.2 仍写"部分" |
| 13 | EF（默认） | §3 | E.3 water_opt −57.79968（仅能量） | A-部分 | 未对照最终几何/梯度范数 |

### 2.2 后处理分支（任一主模式结束后执行）

| 分支 | 文档公式/章节 | 文档中的验证方法 | 判定 | 备注 |
|---|---|---|---|---|
| POLAR | §12 | V3–V5、§19.2（α 9.67480、γ 34.25135、β 20.8707、单位换算） | A | 数值+容差明确 |
| STATIC | §12 | I.8（重定向能量、α 0.8419，修复前后+官方） | A-部分 | 未给 β/γ 全张量 |
| ESP | §13 | E.3 water_esp（主能量全同）；.esp 数值「待逐点比」 | B/部分 | 无 .esp 逐点参考值 |
| PMEP | §13 | H.2「部分：AM1 下 HOF/IP 一致」 | B | 无具体数值表/容差 |
| MULLIK / LOCAL / ENPART / BONDS | §14 | H.2「HOF 一致；分析数值待逐项比」 | B | 仅状态 |
| ALLBO | §14 | H.2「扩展/别名」 | C | 官方不识别，无需对照 |

### 2.3 模式开关与方法（横切分支）

| 模式 | 文档位置 | 文档中的验证方法 | 判定 | 备注 |
|---|---|---|---|---|
| MOZYME（1SCF） | §10 注、§19 | E.3/V25 水 −57.78972 | A | 仅 1SCF |
| MOZYME 几何优化 | 无专章 | E.1 末行「大体系/MOZYME 梯度：待验证」 | C/B | 文档没有优化验证方法 |
| COSMO（1SCF） | §11 | E.3/V23 水 −65.16522 | A | 仅 1SCF |
| COSMO 几何优化 | 无专章 | 无 | C | 文档未涉及 |
| MOZYME + COSMO（linear_cosmo） | I.1 一行（写"未实现"） | 无 | C | 文档写"不可用"，但 2026-10-26 已修复并验证 |
| UHF/ROHF | §10 | E.3 O₂ 7.4358 | A | 仅 O₂ |
| CI/MECI 能量 | §10 | E.3 水 C.I.=2 −58.55115 | A | 仅能量 |
| CI/MECI 分析（getmic/自旋密度/ESR） | I.1「部分」 | 无 | B/C | 部分实现，无验证值 |
| PM7-TS | 方法分支（附录 C/H.1） | E.3 H₂ −83.04626 | A | 仅一个算例 |
| d 元素（Fe/Zn） | §16/C | E.3 FeO/ZnO HOF | A | 仅两种元素 |
| f 元素（La） | §16/C | E.3 LaO HOF | A | 仅一个元素 |
| 6 方法 × s/p（MNDO/AM1/PM3/PM6/PM7/RM1） | 附录 C/D | D.2 37 组 HOF/E_elec/E_nuc 全同 | A | 文档中最完整的矩阵 |
| EXTERNAL / Sparkles | H.1/G.5 提及 | 无 | C | 无验证方法 |
| 周期性 id=1/3（聚合物/晶体） | §15 | E.1「待验证」 | B | 无 tvec deck/参考值 |
| 校正项（色散/氢键/PM6/PM7 专用） | §16 | E.1「隐含验证：部分」 | B | 无逐项（D3/H4/卤键等）对照 |
| 分析输出（DENOUTF/AIGOUT/GRAPH/DIPOLE 等） | H.1 仅列出 | 无 | C | 无输出数值对照 |
| PKA / RAPID / 蛋白-PDB 流程（PDBOUT/RESEQ/RESID/ADD-H/SITE=）/ HTML / FIELD | 无 | 无 | C | 代码有分支，文档完全没有 |

## 3 文档自身声明与统计

- E.1「分支/模式验证状态矩阵」共 23 行：**已验证 10、部分 4、待验证 7、未完成/待修 2**；
- H.2「主要计算模式验证状态」共 24 行：**已验证 9、部分 6、待验证 5、不可用 2、扩展 2**；
- E.5 写"先覆盖 E.1 的 22 类分支"，但 E.1 实为 23 行（统计口径不一致）；
- D.4/E.4/H.3/H.4/I.4 均声明"不能声称全部验证"、"全部验证目前不成立"。

## 4 状态过时与自相矛盾

1. **附录 I 的修复记录没有回写到 E.1/H.2**：
   - E.1 写"FORCE/THERMO 未完成"、H.2 写"不可用"，但 I.7/I.11 已修复（ZPE 10.118、
     熵 45.0697）；
   - E.1 写"STATIC 待修"，但 I.8 已修复（α 0.8419）；
   - E.1 写"L-BFGS-B 待验证"、H.2 写"部分"，但 I.9 已可用（−57.79962 vs −57.79957）；
   - E.1 写"BFGS/DFP 部分"，I.6 已修复（−57.79968 / −57.79966）。
2. **文档早于 2026-10-26 的修复**：MOZYME 几何优化、MOZYME+COSMO、COSMO 几何优化、
   COSMO 优化崩溃修复均未收录；这些分支在 D: 工程已用官方双程序验证（例如：
   水 MOZYME 优化首轮梯度 3.818、HOF −57.79909；水 +COSMO 7.768 / −65.21180；
   苯 +COSMO 20.182 / 20.95883；萘、环戊二烯同样一致）。
3. **审计基线变化**：H.1 记录"305 个关键字"；当前仓库 382 个 `.cpp` 中不同
   `keywrd.find` 字面量已达 **379** 个，文档的关键字审计清单已不是当前状态。
4. **未引用仓库验证资产**：`docs/_verify_runs/` 下已有 11 组 `branches`、23 组
   `options`、42 组 `matrix`、24 组 `regress`（cpp/off 双跑）等实测输出，以及
   `run_gate.py`、`run_cosmo_moz.py`、`run_moz_opt_probe.py` 三个回归脚本；
   文档只使用其中一部分结果，且没有给出目录/脚本索引。

## 5 建议补齐（按优先级）

1. **新增附录 J「逐分支验证方案表」**，列：分支 / 算例 deck / 关键字 / 对照量 /
   参考值+容差 / 脚本或输出目录 / 状态。以 `docs/_verify_runs/` 为执行载体。
2. **回写并统一状态**：E.1、H.2、G.2、E.3 与 I.6–I.11 及 2026-10-26 修复对齐；
   明确"已验证"分级（仅能量 / 能量+几何 / 全输出）。
3. **补齐 B 类（无具体方法）** 的最小双程序算例：
   TS/SADDLE/NLLSQ/SIGMA（HCN/HNC）、PATHS/PATHK/GRID（2–3 点扫描）、
   LOCATE/REFINE-TS（小蛋白）、周期性（tvec polymer/crystal）、DRC 轨迹逐点、
   PMEP 数值表、ESP .esp 逐点、分析输出（MULLIK/LOCAL/ENPART/BONDS）逐项。
4. **补齐 C 类**：PKA、RAPID、EXTERNAL、Sparkles、蛋白/PDB 流程、HTML、FIELD、
   DIPOLE 输出、MOZYME 优化、COSMO 优化、MOZYME+COSMO、校正项逐项。
5. 把 `run_gate.py`（14 项）、`run_cosmo_moz.py`（12 项）、`run_moz_opt_probe.py`
   （15 项）作为附录 J 的执行入口，并在文档中给出当前通过的基线值
   （THERMO ZPE 10.007/官方 10.003、FORCE ZPE 10.120/官方 10.119、COSMO 套件
   10/12 完全一致等）。

## 6 文档更新记录（2026-10-26）

本审计已落实为 `FORMULAS_AND_VERIFICATION.docx` 的 **附录 J**（J.1–J.8，5 张表），
并同步修改原文档：

- **回写 E.1**（附录 E 状态矩阵，15 行状态/对照量更新）与 **H.2**（附录 H 状态表，
  16 行更新）：按 I.6–I.11 与 2026-10-26 的 MOZYME/COSMO 修复同步，消除
  E.1/H.2 与 I.6–I.11 的自相矛盾；
- **验证索引新增 V47–V50**：MOZYME 几何优化、COSMO 几何优化、MOZYME+COSMO、
  分支覆盖审计；
- E.1/H.2 表后各加一条 **2026-10-26 修订说明**，指向附录 J；
- 原有 72 个 Word 原生公式（OMML）与排版保持不变（更新前后 OMML 计数 72/72，
  表格 37→42，仅新增附录 J 的 5 张表）；
- 修改前的文件备份：工作区 `port_chklew/FORMULAS_AND_VERIFICATION.backup-20261026.docx`。

后续每次修复分支后，按 J.7 的 P0–P3 顺序更新 E.1/H.2/J 与验证索引。

## 7 按“读者可独立复算”目的的复核（2026-10-26，附录 K）

文档用途进一步明确为：读者根据文档的公式与参数，自己手算或编程算出程序结果，
用于验证并判断实现是否合格。按此目的复核（判定五要素：输入条件 / 完整公式 /
参数 / 中间参考值 / 最终参考值+容差），结论是**目前只对少数分支成立**：

| 可复算级别 | 分支/内容 |
|---|---|
| L1 可手算 | §20 的 10 个案例（BFGS/DFP、P-RFO、Verlet、CPHF、CI 2×2、配分函数、色散/COSMO、ZPE 等）；附录 A.5 的 H₂ SCF 骨架 |
| L2 自编程复算 | 附录 A/B/F 的 H₂/PM7 全链（含中间值）；附录 D 的 37 组 1SCF 最终值；d/f 已测组合；FLEPO 单步等 |
| L3 需程序中间量 | EF、FORCE/THERMO、POLAR/STATIC、COSMO、CI/MECI、UHF、PM7-TS、L-BFGS-B |
| L4 不可独立复算 | TS/SADDLE/NLLSQ/SIGMA、LOCATE/REFINE-TS、PATHS/PATHK/GRID、DRC/MOLDYN、周期性、ESP/PMEP 细节、分析输出、MOZYME、MOZYME/COSMO 优化、MOZYME+COSMO、PKA/RAPID/蛋白-PDB/HTML/FIELD/EXTERNAL/Sparkles |

关键缺口（详见文档附录 K.3）：

1. 中间参考值仅 H₂ 全链 + §20 手算案例 + I.6–I.11 修复记录；其余分支无逐步骤
   H/W/F/P/ε/E_elec/E_nuc 可对照；
2. 附录 F 只给 s-s 重叠闭式，s/p/d 重叠与双电子/γ/ecore 多为形式或源码位置；
3. 多数分支没有可复制的 deck（关键字+几何+电荷/自旋）；
4. 参数表缺 COSMO 半径/方向、色散/氢键、PM7-TS、MOZYME 阈值等；
5. 默认/收敛常量（scfcrt、GNORM、DMAX、trust radius、DIATOM/DIIS 等）未集中列出；
6. 参考输出+容差只覆盖部分分支，L4 分支无法判定“合格”；
7. MOZYME/linear_cosmo/COSMO-MOZYME 只有 2026-10-26 的最终数值，无公式与中间量。

已落实的文档更新（2026-10-26，第二次更新）：

- 文档开头新增“按读者独立复算目的”的说明，指向附录 K；
- 新增 **附录 K**（K.1–K.7，4 张表）：五要素/四级判定、逐分支可复算级别与缺失要素、
  关键缺口、补齐方案、结论、**水 PM7 1SCF 的复算检查点示例**（HOF/E_elec/E_nuc/
  本征值/净电荷/布居/偶极 + 容差）、每分支检查点表字段模板；
- 附录 J.8 增加指向 K 的说明；
- 公式（OMML 72 个）保持不变；表格 42→46（仅新增 K 的 4 张表）；备份文件：
  `port_chklew/FORMULAS_AND_VERIFICATION.backup2-20261026.docx`。

下一步（使文档达到其目的）：按 K.7 模板，用 `docs/_verify_runs` 的双程序输出
脚本化生成每个分支的“复算检查点表”（deck + 中间量 + 最终值 + 容差），并补齐
K.3 的 7 类缺口；优先从 L2/L3 分支开始（1SCF/EF/FORCE/POLAR/COSMO/UHF/CI）。

## 8 第一批复算检查点表（2026-10-26，已写入文档附录 K.8）

从 `docs/_verify_runs` 的官方输出（ref）与 C++ 输出（cpp）中提取并写入
`FORMULAS_AND_VERIFICATION.docx` 的 **K.8**（4 张表）：

1. **K.8.1 能量与电子结构**（14 行）：1SCF/EF/FLEPO/DFP/L-BFGS/COSMO/UHF/CI/
   MOZYME/PM7-TS/ZnO(PM6/PM7)/LaO/H₂ 的 ΔHf、TOTAL、E_elec、E_nuc、IP、
   HOMO/LUMO、⟨S²⟩、STATE 等，并注明 C++ 对照；
2. **K.8.2 优化器迭代轨迹**（4 行）：EF/FLEPO/DFP/L-BFGS 的逐轮（梯度/ΔHf）、
   终点 GNORM 与 O–H/H–H 键长检查点；
3. **K.8.3 FORCE/THERMO**（10 行）：FORCE ZPE 10.119/10.120、三频率、
   笛卡尔力常数；THERMO ZPE 10.003/10.007、298K 焓/Cp/熵、ROT/TRA/VIB 分项；
4. **K.8.4 响应与分析输出**（10 行）：STATIC α、MULLIK、LOCAL、ENPART、BONDS、
   ESP、PMEP、IRC、DRC、0SCF 的参考值与 C++ 状态。

同时记录本次提取中发现的**已知偏差/缺口**（已写入 K.8.3/K.8.4/K.8.5）：

- EF/FLEPO/DFP/L-BFGS 的能量与 ΔHf 对到打印精度，但迭代轨迹/终点 GNORM 不同
  （EF 第 2 轮 C++ 1.726 vs 官方 1.741；L-BFGS 阈值 0.5 vs 1.0 等）；
- C++ FORCE 的频率打印顺序/标签异常（ZPE 正确：10.120 vs 10.119），
  FORCE CONSTANT 块未输出；
- C++ THERMO 的 ZPE/焓/Cp/熵有 10⁻³~10⁻¹ 差，频率集与 FORCE 不同（I.11 已记）；
- C++ 的 LOCAL/BONDS/ESP/PMEP 打印尚未实现（不影响能量）；
- IRC 初速轴序不同、逐点能量差 10⁻³~10⁻²；DRC 轨迹点数/长度不同；
- 0SCF 仅有结构输出，无能量参考（H.2 记为“输出结构”）。

文档当前结构：正文/附录 A–J + K.1–K.7（可复算性评估）+ **K.8（第一批检查点表）**；
表格 46→50，公式（OMML）仍为 72 个不变；备份文件：
`port_chklew/FORMULAS_AND_VERIFICATION.backup3-20261026.docx`。

下一批（待做）：把 K.8 中“待补/部分”的项目补成可复算数据（H/W/F/P/ε 中间量
dump、FORCE 频率与力常数打印、LOCAL/BONDS/ESP/PMEP 输出、IRC/DRC 轨迹对齐），
并按 K.8 的格式扩展 TS/PATHS/GRID/周期性等 L4 分支。
