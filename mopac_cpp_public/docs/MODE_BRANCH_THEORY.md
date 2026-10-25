# MOPAC-CPP 各模式分支的数学理论

> 本文按 **run_mopac 的计算分派分支** 逐一给出该模式的数学问题、控制方程、
> 迭代/求解步骤、收敛判据与输出量。所有公式与 C++ 源码位置对应；
> 底层 SCF/NDDO 积分与 MOZYME 伪对角化的推导细节见
> `COMPUTATION_FORMULAS.md`，本文引用其记号，不重复推导。

---

## 0. 记号与公共基础

| 记号 | 含义 |
|---|---|
| **R** | 3N 维核坐标（内部或笛卡尔） |
| P | 密度矩阵（打包下三角，`j(j-1)/2+i`） |
| H, F | 单电子矩阵、Fock 矩阵 |
| W | NDDO 双电子积分 |
| C, ε | MO 系数、轨道能 |
| E_elec, E_nuc, ΔHf | 电子能、核-核排斥能、生成热 |
| g = ∂ΔHf/∂**R** | 解析/数值梯度 |
| **H** = ∂²ΔHf/∂**R**² | Hessian（质量加权或普通） |
| gnorm = |g| | 梯度范数 |
| fpc₉ = 23.060529 | eV → kcal/mol |

所有分支共用的单点能泛函：

```
E(R) = min_{P} { ½ Σ_μν P_μν (H_μν(R) + F_μν(P,R)) + E_nuc(R) } + E_corr(R)
ΔHf(R) = fpc₉ · (E_elec + E_nuc) + atheat + E_disp + E_hb + E_solv + …
```

其中 `E_corr` 含色散、氢键、PM6/PM7 专用校正、COSMO 溶剂化、压力项等
（见 `COMPUTATION_FORMULAS.md` §7）。SCF 内层在 **每个模式** 中都是
“Fock→对角化→密度→收敛”的同一套迭代（非 MOZYME）或
“buildf→diagg1/diagg2→density”的 MOZYME 迭代。

---

## 1. 分支总表（`run_mopac.cpp` 分派顺序）

| 优先级 | 关键字/条件 | 模式 | 入口 | 目标 |
|---|---|---|---|---|
| 1 | `SADDLE` | 鞍点搜索 | `react1()` | 沿两构型连线找能量极大 |
| 2 | `STEP1`（`GRID` 关键字设置） | 二维网格 | `grid()` | 势能面 (E vs 两个内坐标) |
| 3 | `latom != 0` | 路径跟随 | `paths()`/`pathk()` | 固定反应坐标扫描 |
| 4 | `FORCE` / `IRC=` / `THERM` / `DFORCE` | 力常数 | `force()` | Hessian、频率、热力学 |
| 5 | `DRC` / `IRC` | 动力学 | `drc()` | 经典轨迹 / IRC |
| 6 | `MOLDYN` | 分子动力学 | `run_md_velocity_verlet()` | NVT/NVE 轨迹 |
| 7 | `NLLSQ` | 最小二乘过渡态精修 | `nllsq()` | 二次型拟合 |
| 8 | `SIGMA` | σ 键约束精修 | `sigma()` | 约束优化 |
| 9 | `1SCF` 或 `nvar==0` | 单点 | `compfg(x,true,…)` | E(R) 一次 SCF |
| 10 | `LOCATE-TS` / `REFINE-TS` | 过渡态 | `Locate_TS_for_Proteins()` / `Refine_TS_for_Proteins()` | 蛋白/大体系 TS |
| 11 | `DFP`/`FLEPO`/`BFGS`/`nvar==1` | BFGS/DFP 优化 | `flepo()` | min E(R) |
| 12 | `TS` | EF 过渡态 | `ef()` | 鞍点（max 1 模） |
| 13 | `LBFGS` / `GEO_REF` / `nvar>1000` | L-BFGS-B | `lbfgs()` | min E(R) |
| 14 | 默认 | EF 优化 | `ef()` | min E(R) |
| 后处理 | `POLAR` / `STATIC` / `PMEP` / `ESP` / `MULLIK` / `SYB` / `LOCAL` / `ENPART` / `BONDS` | 响应/电荷/分析 | 各自模块 | 见 §12–§16 |

修饰性分支（与上述任一模式正交）：
`MOZYME`（LMO SCF）、`UHF`/`CI`/`ROOT`（电子结构）、
`id=1/3`（一维/三维周期）、`EPS`/`COSMO`（溶剂化）、
`0SCF`（只做 setup）、`DISP`/`PM6-D`/`PM7`（校正项）。

---

## 2. 0SCF：只做数据准备（无 SCF）

- 触发：`0SCF`。
- 数学内容：无迭代；只解析关键字/几何、调用 `moldat` 建立电子数与原子参数、
  `hcore` 生成 H/W/E_nuc，然后 `writmo` 输出初始结构。
- 输出：`E_nuc`（核-核排斥）、初始几何、经验式、点群。

---

## 3. 1SCF：单点能

- 触发：`1SCF` 或优化变量数 `nvar==0`。
- 数学问题：

```
E(R₀) = min_P ½ Σ_μν P_μν (H_μν + F_μν(P)) + E_nuc
```

- 求解：固定 R₀，跑一次完整 SCF（`iter()`/`iter_for_MOZYME()`），
  只做能量与可选梯度，不移动核。
- 收敛：SCF 密度判据 `pl < scfcrt`（POLAR 时 `scfcrt ≤ 1e-6`，
  见 `iter.cpp`）。
- 输出：ΔHf、轨道能、偶极、电荷、选定的打印块。
- 位置：`run_mopac.cpp:710-721`、`compfg.cpp`。

---

## 4. EF：本征向量跟随优化（默认）

- 触发：默认分支、`TS`、`EF`。
- 数学问题：在 Hessian 本征基下分别控制“极小化/极大化”方向。

设 `H v_i = λ_i v_i`，梯度投影 `g_i = v_iᵀ g`。对极小点，
位移满足：

```
Δx = - Σ_i  (g_i / (λ_i - λ)) v_i           （P-RFO）
```

λ 是移位参数（trust-radius 拉格朗日乘子），由 `|Δx| = R_trust` 确定。
对过渡态（TS），沿跟踪模式 m 用“最大化”移位 λ_p（使 λ_m - λ_p > 0），
其余模式用“最小化”移位 λ_n（λ_i - λ_n > 0）：

```
Δx = - (g_m/(λ_m - λ_p)) v_m
     - Σ_{i≠m} (g_i/(λ_i - λ_n)) v_i
```

λ_p、λ_n 各自通过一维求根使两个子空间位移满足同一 trust radius。
`formd()` 实际实现三种步型：
- **NR**（纯 Newton-Raphson，`Δx = -H⁻¹g`，仅作参考/调试）；
- **P-RFO**（上述分区移位，极小点用同一 λ）；
- **QA**（对 λ 做二次近似，等价于 trust-region 的近似解）。

极小点用单一 λ 使 |Δx|=R_trust；TS 用 λ_p、λ_n 两个移位分别满足
投影长度约束。

- 迭代步骤（`ef()` 主循环）：
  1. `compfg` 得 E、g；
  2. Hessian 本征分解（`eigenvectors_LAPACK`/`rsp`）；
  3. `overlp()` 把新本征向量与旧模式做重叠匹配（mode following，
     防本征值交叉导致跳变）；
  4. `prjfc()` 投影掉平动/转动分量；
  5. `formd()` 求位移 Δx（trust radius 控制）；
  6. 更新坐标，必要时 `updhes()` 做 Powell/BFGS Hessian 更新；
  7. 收敛判断，未收敛回 1。
- 收敛判据：|g|（GNORM）、|ΔE|、|Δx|、模式跟踪状态、Herbert 检验。
- 初始 Hessian：`gethes()` 对梯度做数值差分得到。
- 输出：优化后的 ΔHf、几何、频率（若同时 `FORCE`）。
- 位置：`ef.cpp`、`ef_C.cpp`；分派 `run_mopac.cpp:740-752`。

---

## 5. FLEPO：BFGS/DFP 变尺度优化

- 触发：`DFP`/`FLEPO`/`BFGS` 或 `nvar==1`。
- 数学问题：逆 Hessian 近似 **B** ≈ **H**⁻¹，方向 `p = -B g`。

BFGS 更新（默认）：

```
s = x_{k+1} - x_k ,   y = g_{k+1} - g_k ,   ρ = 1/(yᵀ s)
B_{k+1} = (I - ρ s yᵀ) B_k (I - ρ y sᵀ) + ρ s sᵀ
```

DFP 变体（关键字 `DFP`）：

```
B_{k+1} = B_k + (s sᵀ)/(sᵀ y) - (B_k y yᵀ B_k)/(yᵀ B_k y)
```

- 线搜索 `linmin()`：沿 p 最小化 φ(α)=ΔHf(x+αp)，用二次/三次插值 +
  梯度内积判据，得到可接受 α；失败按 `flepo` 的回退逻辑缩步。
- 收敛：Peters 检验（连续能量变化同号且小）、|g|、|Δx|、ΔE；
  线搜索失败两次会给出警告信息。
- 重启：`dfpsav()` 周期保存 Hessian 与几何。
- 位置：`flepo.cpp`、`linmin`、`dfpsav`。

---

## 6. L-BFGS-B：有限内存拟牛顿（大体系）

- 触发：`LBFGS`、`GEO_REF`（id=0）、`nvar>1000` 且未指定 EF。
- 数学问题：min E(x)，x 满足关键字给定的上下界（投影梯度）。

两循环递归求方向（m 个存储对 `{s_i, y_i}`，`ρ_i=1/(y_iᵀ s_i)`）：

```
q = g
for i = k-1 … k-m:  α_i = ρ_i (s_iᵀ q);  q -= α_i y_i
r = H₀ q ,  H₀ = γ I ,  γ = (s_{k-1}ᵀ y_{k-1})/(y_{k-1}ᵀ y_{k-1})
for i = k-m … k-1:  β = ρ_i (y_iᵀ r);  r += s_i (α_i - β)
p = -r
```

- 投影/活动集：把 x 限制在 `[lo, hi]`，按广义 Cauchy 点与子空间最小化求步长；
  实现见 `lbfgs.cpp` 的 `mainlb/subsmem/matupd` 等。
- 收敛：投影梯度范数 ≤ `pgtol`、|ΔE| ≤ `factr·epsmch`、最大迭代。
- 输出：优化几何、ΔHf、可选 Hessian 近似。
- 位置：`lbfgs.cpp`、`lbfgs_line_search.cpp`。

---

## 7. 过渡态各分支

### 7.1 `TS`（EF 鞍点优化）
与 §4 相同，但跟踪**最低非零本征值**模式并沿其最大化：λ_p 取在
λ_m 之上的移位，其余模式用 λ_n 极小化（P-RFO）。输出虚频/鞍点几何
（配合 `FORCE` 验证 1 个虚频）。

### 7.2 `NLLSQ`（最小二乘精修）
- 数学问题：把梯度分量作为残差 r(x)=g(x)，最小化 ½|g|²（等价于
  在 TS 附近精修），用 Gauss-Newton 型步：
  `Δx = -(JᵀJ)⁻¹ Jᵀ g`，J=∂g/∂x（Hessian 近似）。
- 收敛：|g| 下降、能量变化、最大迭代。

### 7.3 `SIGMA`（σ 约束精修）
- 在若干键长/键角上加强约束，用惩罚/约束最小二乘修正梯度方向，
  保持指定 σ 骨架不变（见 `sigma` 实现与关键字解析）。

### 7.4 `SADDLE`（react1 双构型搜索）
- 维护两套几何 A、B；各自 SCF 得 ΔHf_A、ΔHf_B；
- 将较低能量体系沿 A→B 连线方向移动一步：
  `x_high ← x_high + α (x_low - x_high)`（α 由关键字控制），
  重新优化/SCF，直到能量差或距离收敛，得到鞍点近似。
- 位置：`react1.cpp`。

### 7.5 `LOCATE-TS` / `REFINE-TS`（蛋白/大体系）
- LOCATE-TS：在给定反应坐标约束下先做一次过渡态定位（RAPID 模式下用
  Hessian 最小模式方向爬升），再做 EF/NLLSQ 精修；
- REFINE-TS：用已有近似 TS 做 P-RFO 精修，限制重原子移动（大体系）。
- 位置：`Locate_TS_for_Proteins.cpp`、`Refine_TS_for_Proteins.cpp`。

---

## 8. 反应路径与扫描

### 8.1 `PATHS`（固定反应坐标扫描，`paths.cpp`）
- 数学问题：对反应坐标 r(x) 取离散值 `r_1<…<r_K`，每个点解
  `min E(x) s.t. r(x)=r_k`（其余变量优化，默认 EF）。
- 实现：逐步改变反应坐标值 → 热启动 EF/DFP → `writmo` 输出该点；
  可生成轨迹/PDB。
- 输出：E(r) 曲线、每点几何。

### 8.2 `PATHK`（等步长路径跟随，`pathk.cpp`）
- 沿路径以固定步长（关键字 `STEP`、`POINTS`）推进：在每步局部坐标系中
  把切向固定、优化垂直方向；用于 IRC 近似或反应路径。

### 8.3 `GRID`（二维势能面，`grid.cpp`）
- 两个内坐标 (v₁,v₂) 与范围/步长；网格点 `(v₁_i, v₂_j)`：
  `E_ij = min_{x⊥} E(x)`（或 fixed-point SCF）。
- 蛇形四遍扫描（左→右、右→左、上→下、下→上），相邻点热启动，
  每点记录最低能量几何；输出 `.grid`/`.arc` 表。

---

## 9. 力常数、振动频率与热力学

### 9.1 FORCE/DFORCE：Hessian 与频率（`force.cpp`）
- 力常数（Hessian）由**解析梯度**做数值差分得到（`compfg` 在位移构型上
  重新算梯度）：
  `H_ij ≈ [g_i(x+h e_j) - g_i(x-h e_j)]/(2h)`；
  `DFORCE` 是调试/直接力常数路径，与 `FORCE` 共用同一套位移-梯度框架。
- 质量加权 Hessian：`H̃ = M^{-1/2} H M^{-1/2}`，`M_ii = m_atom`。
- 本征方程：`H̃ φ_k = λ_k φ_k`。
- 简谐频率：`ν_k = (1/2πc) √(λ_k)`（λ<0 → 虚频，TS 应恰有 1 个）。
- 零点能：`ZPE = ½ Σ_k h ν_k`。
- 同位素替换：改变 M 后重解特征问题。
- 输出：频率表、约化质量、力常数、IR 强度（若有关键字）、ZPE。
- 位置：`force.cpp`。

### 9.2 THERMO：统计热力学（`thermo.cpp`）
- 平动配分函数：`Q_tr = (2π m kT/h²)^{3/2} V`
- 转动（线性）：`Q_rot = 8π² I kT/(σ h²)`；
  非线性：`Q_rot = √π/σ · (8π² kT/h²)^{3/2} √(I_A I_B I_C)`
- 振动：`Q_vib = Π_k (1 - e^{-hν_k/kT})^{-1}`
- 电子：`Q_elec = g_elec`（简并度）
- 总配分函数：`Q = Q_tr Q_rot Q_vib Q_elec`
- 内能/焓：`U = kT² ∂lnQ/∂T`，`H = U + kT`（理想气体）
- 熵：`S = k lnQ + kT ∂lnQ/∂T`
- 热容：`C_p = (∂H/∂T)_p`，`C_v = (∂U/∂T)_V`
- 生成热温度修正：由 298 K 的 ΔHf 与原子化能/热容积分给出
  `ΔHf(T) = ΔHf(298) + ∫_{298}^{T} ΔC_p dT`；
  MOPAC 输出 H、Cp、S、G 表（`THERMO` 控制温度范围）。
- 位置：`thermo.cpp`。

---

## 10. 动力学与反应坐标

### 10.1 DRC / IRC（`drc.cpp`）
- 经典运动方程：`M ẍ = -∇E(x)`（内部/笛卡尔坐标，含 G 矩阵变换）。
- 速度 Verlet：
  `x_{n+1} = x_n + v_n Δt + ½ a_n Δt²`
  `v_{n+1} = v_n + ½ (a_n + a_{n+1}) Δt`，`a_n = -M⁻¹ ∇E(x_n)`
- 初始化：由 `startv`/`startk` 给定速度并按质心平动/转动归零。
- 阻尼/温控：按关键字对动能乘阻尼因子（DRC 常用）。
- 终止：步数/时间上限、能量阈值、超时写 `.res` 重启。
- 输出：轨迹（`.arc/.xyz`）、每步 E、梯度范数、动能。
- 位置：`drc.cpp`。

### 10.2 MOLDYN（速度 Verlet MD，`run_md_velocity_verlet`）
- 触发：`MOLDYN` 关键字，或注释行含 `Mixture:`（`run_mopac.cpp:691-700`）。
- 与 10.1 同一积分器，另加恒温器（Berendsen/Nose 由关键字选择）：
  `λ = √(1 + (Δt/τ)(T₀/T - 1))`，`v ← λ v`。
- 输出：时间序列能量/温度/压力；轨迹。

---

## 11. 电子结构分支

### 11.1 RHF/UHF/ROHF
- RHF：`F = H + G(2P)`，`P = Σ_occ C Cᵀ`。
- UHF：`F^α = H + G(P^α)`, `F^β = H + G(P^β)`；
  `P = P^α + P^β`，自旋密度 `P^s = P^α - P^β`；
  `E = ½ Σ (P^α:F^α + P^β:F^β)`；
  `S_Z = (N_α - N_β)/2`，`⟨S²⟩ = S_Z² + Σ |P^s|²` 修正。
- 开壳分数占据：占据数 `n_i = 1 ± fract`，Fock 取 `F = H + G(P)`，
  密度按分数组合。

### 11.2 MECI / CI（`meci.cpp`）
- 活性空间：选 `nmos` 个轨道与 `ielec` 个电子，枚举 α/β 微态（`perm`）。
- CI 矩阵：`H_{IJ} = Σ_{ij} h_ij γ^I_{ij} + ½ Σ_{ijkl} (ij|kl) Γ^I_{ijkl}`
  （γ、Γ 为单/双占据数矩阵，由微态生成）。
- 对角化：`H C = E C`（小空间直接 LAPACK，大空间 Davidson 迭代）。
- 取目标根 `E_root`，CI 能量修正加回总能量；输出 CI 系数、自旋密度、ESR 量。
- 位置：`meci.cpp`、`mecid.cpp`、`mecih.cpp`。

### 11.3 MOZYME（与上述任一模式正交）
- 单点/优化/频率都可用 `MOZYME`：LMO SCF（`diagg1/diagg2` 伪对角化，
  见 `COMPUTATION_FORMULAS.md` §3），密度 `P = Σ_occ C Cᵀ`，
  能量仍按 §0；梯度用 MOZYME 专用解析式（`prtgra` 输出）。
- 触发：`MOZYME` 关键字；`iter_for_MOZYME` 替代 `iter`。

---

## 12. 溶剂化 / COSMO（`EPS=`/`COSMO`）

- 介电因子：`f(ε) = (ε-1)/(ε+½)`。
- 腔体：SAS 面元面积 `A_i`、体积 `V`；`B` 矩阵（面元间 `1/r` + 对角修正）。
- 屏蔽电荷：`q = -f(ε) · A⁻¹ V`，
  `V_i = Σ_A Z_A/|r_i - R_A| - Σ_μν P_μν (μν|i)`。
- 介电能：`E_diel = ½ Σ_i q_i V_i`。
- 反馈：`F ← F + B·q` 项（`addfck/addnuc`），每轮 SCF 更新。
- 总生成热：`ΔHf += fpc₉ · E_diel`。
- 面积/体积输出：`area`、`cosvol`。
- 位置：`cosmo.cpp`。

---

## 13. 极化率/超极化率

### 13.1 POLAR：TDHF（`polar.cpp`）
- CPHF（含频 ω）：
  `(ε_a - ε_i - ω) u_ai + Σ_bj [2(ai|bj) - (aj|bi)] u_bj = -d_ai^μ`
- α：`α_ij = -Σ_ai d_ai^i u_ai^j`；各向同性 `α = (α_xx+α_yy+α_zz)/3`。
- β（SHG）：二阶响应，`betaf` 迭代 / `nonbet` 非迭代；
  分量平均 `bav = (1/5)Σ`，`bvec = |bav|`。
- γ（THG）：`ngamtg` 9 分量、`gav/5` 平均。
- 迭代：`ATOL/BTOL` 判据，`MAXITA/MAXITU` 上限；频率 `E=(…)`。
- 单位：a.u. ↔ Å³（`a0³`）、ESU（`au_to_esu`）。

### 13.2 STATIC：有限场（`static_polarizability.cpp`）
- 在 ±x,±y,±z 方向加场 F，展开：
  `E(F) = E₀ - Σ_i μ_i F_i - ½ Σ_ij α_ij F_iF_j
        - (1/6) Σ_ijk β_ijk F_iF_jF_k - (1/24) Σ_ijkl γ_ijkl F_iF_jF_kF_l`
- 由差分得到 μ、α、β、γ；张量用 Jacobi 对角化，输出各向同性值。
- 共 6 方向 × 若干场强，每次调用 `compfg`（36 次 SCF 的由来）。

---

## 14. 静电势与电荷

### 14.1 ESP（`esp.cpp`）
- 静电势：`V(r) = Σ_A Z_A/|r-R_A| - Σ_μν P_μν ∫ φ_μ φ_ν/|r-r'| dr'`
  （AO 积分，`elesn`）。
- 采样：Williams 网格（`pdgrid`）、Connolly 表面（`surfac`）等。
- 拟合：min Σ_p [V_p - Σ_A q_A/|r_p-R_A|]²，
  约束 Σ q_A = 总电荷（可选），最小二乘（`espfit`）。
- 输出：`.esp/.esr` 电荷与势。

### 14.2 PMEP（`pmep.cpp`）
- 同 ESP 的 MEP（AM1 参数化），`mepchg` 拟合、`mepmap` 等值线、
  `meprot` 选择映射平面。

---

## 15. 布居/成键分析

- MULLIK：`P_A = Σ_{μ∈A} (P S)_μμ`，净电荷 `q_A = Z_A - P_A`。
- SYB（`mpcsyb`）：Mulliken/SYBYL 格式输出。
- BONDS/ALLBO：键级由 `P S` 与原子对贡献给出。
- LOCAL：局域化轨道（`local`/`local_for_MOZYME`，最大化 Σ|ψ|⁴）。
- ENPART：能量分解（单中心/双中心项）。
- 位置：`mullik.cpp`、`mpcsyb.cpp`、`bonds.cpp`、`local.cpp`、`enpart.cpp`。

---

## 16. 周期性（id=1 聚合物 / id=3 晶体）

- 晶格矢量 `tvec`（3×3）；原胞内坐标 `x`；平移向量 `T` 求和。
- 晶格能量：`E_lattice = ½ Σ_{T} Σ_{A,B} E_AB(R_A - R_B - T)`
  （`hcore` 的 `solrot` 负责周期求和/截断）。
- 压力/应力：`P = -(1/V) ∂E/∂V`，
  `hpress = -P·V`（id=3）或 `-P·|tvec|`（id=1）。
- 晶胞优化：对 `tvec` 的 3×3 分量与内坐标一起优化（`wrtkey` 的 cell 变量）；
  `write_cell/write_pressure` 输出。
- 密度：`ρ = M·10²⁴/(N_A·V)`（`write_cell` 公式）。

---

## 17. 校正项（所有模式共用）

| 项 | 公式 |
|---|---|
| 色散（PM6-D/PM7） | `E_disp = Σ_{A<B} f_damp(R) [C₆/R⁶ + C₈/R⁸ + C₁₀/R¹⁰]` |
| 氢键 | `E_hb = Σ htype·sin²θ`（`addhb`/H_bond 模块） |
| PM6 >N− | `nsp2_correction()`，加进 atheat |
| PM7 Si-O-H | `Si_O_H_Correction()` |
| 二面角 | `sum_dihed = Σ htype·sin²θ` |
| 压力 | `±p·V` |
| 溶剂化 | COSMO `E_diel` |
| 外场 | `FIELD` → `fldon` 加进 H/Fock |

---

## 18. 各模式的收敛/终止判据

| 模式 | 判据 |
|---|---|
| SCF | `pl < scfcrt`（POLAR 时 ≤1e-6），ΔE 小 |
| EF | |g|、|ΔE|、|Δx|、模式跟踪、Herbert |
| BFGS/DFP | Peters 检验、|g|、|Δx|、ΔE，线搜索成功 |
| L-BFGS-B | 投影梯度 ≤ pgtol、|ΔE| ≤ factr·epsmch |
| SADDLE | 两构型能量差/距离收敛 |
| PATHS/GRID | 指定步数/范围走完 |
| DRC/MD | 步数/时间上限、能量阈值 |
| FORCE/THERMO | 线性方程组规模固定，一次求解 |
| POLAR | |Δu| < ATOL/BTOL，迭代 ≤ MAXIT |
| CI | Davidson 残差/根数 |
| 几何优化 | 由优化器判据决定（上级 §4–§6） |

---

## 19. 单位与常量

- 能量：内部 eV；生成热 kcal/mol（`fpc₉=23.060529`）。
- 长度：Å（`a0=0.5291772083`）。
- 频率：cm⁻¹；`ν̃ = (1/2πc)√λ`。
- 极化率：a.u.（1 a.u. = a0³ Å³）、ESU（`au_to_esu=8.639418e-33`）。
- 温度：K；`k_B`、`h`、`c` 见 `funcon_C`。

---

## 20. 文档与源码对照

| 主题 | 文档 | 源码 |
|---|---|---|
| NDDO/SCF 基础 | `COMPUTATION_FORMULAS.md` §1-2 | `hcore.cpp`, `fock2.cpp`, `iter.cpp` |
| MOZYME | `COMPUTATION_FORMULAS.md` §3 | `iter_for_MOZYME.cpp`, `diagg.cpp` |
| EF/BFGS/L-BFGS | 本文 §4-§6 | `ef.cpp`, `flepo.cpp`, `lbfgs.cpp` |
| TS/路径/网格 | 本文 §7-§8 | `react1.cpp`, `Locate_TS…`, `paths.cpp`, `grid.cpp` |
| 频率/热力学 | 本文 §9 | `force.cpp`, `thermo.cpp` |
| 动力学 | 本文 §10 | `drc.cpp`, `run_md_velocity_verlet` |
| CI/UHF | 本文 §11 | `meci.cpp`, `iter.cpp` |
| COSMO | 本文 §12 | `cosmo.cpp` |
| POLAR/STATIC | 本文 §13 | `polar.cpp`, `static_polarizability.cpp` |
| ESP/PMEP | 本文 §14 | `esp.cpp`, `pmep.cpp` |
| 布居分析 | 本文 §15 | `mullik.cpp`, `bonds.cpp` |
| 周期性 | 本文 §16 | `hcore.cpp`, `write_cell.cpp` |
| 校正项 | 本文 §17 | `compfg.cpp`, 各 correction 模块 |

> 说明：本文给出各模式的“数学骨架”，具体数值实现细节（循环顺序、
> 打包约定、判据常量）以对应 `.cpp` 的中文导读与源码为准；
> 位级可复现性边界见 `BITWISE_FEASIBILITY.md`。

---

## 21. 公式验证（2026-10-25 增补）

完整版 Word 文档 `MODE_BRANCH_THEORY.docx`（13 页、54 个公式对象）新增：

- 全局符号表（符号 → 含义 → 单位/取值）；
- 每个模式的 **输入 / 公式 / 逐步迭代 / 输出 / 收敛判据 / 源码**；
- 第 19 章公式验证：用真实输出核对
  `TOTAL=ELEC+CORE`、`ΔHf−fpc₉(E_elec+E_nuc)=atheat_eff`、`α_iso=mean(diag)`、
  `γ_avg=Σγ/5`、`|b|=√Σb²`、`ZPE=½Σhcν`、`α(Å³)=α(a.u.)a₀³`；
  并用随机矩阵/ODE 独立校验 BFGS/DFP 割线条件、P-RFO→NR 极限、
  速度 Verlet 能量守恒、CPHF 标量解、配分函数量级。
- 验证中发现并修复 `force.cpp` 的 `xyzint` 指针偏移错误；记录 `fmat()` 占位
  （FORCE 频率为 0）与 STATIC SCF 不收敛两个剩余缺口。

验证结果详见 `TRANSLATION_GAPS.md` §57。

---

## 22. 手工验算（2026-10-25 增补）

Word 文档新增第 20 章：基础参数表（PM7 的 eheat/eisol、水的 E_elec/E_nuc、
常数）、10 个可手工复算案例（生成热、POLAR、ZPE、BFGS、P-RFO、Verlet、
CPHF、2×2 CI、配分函数、色散/COSMO）及验算清单（应得值/容差）。
详见 `TRANSLATION_GAPS.md` §59。

---

## 23. 手工验算的可复算范围（2026-10-25 增补）

第 20 章新增 20.0 节：把手工验算分为 L1（完全手算，案例 4–10）、
L2（公式链完整但需程序中间量，案例 1–3）、L3（从几何到 SCF 不可手算）；
并新增 20.13 的 2×2 SCF 教学演示，展示“积分→Fock→对角化→密度→判据”
的完整循环。详见 `TRANSLATION_GAPS.md` §60。

---

## 24. 20.13 完整迭代表（2026-10-25 增补）

20.13 现给出从 H、F=H+0.5P、P₀=0 开始的完整两轮 SCF 迭代、每轮 F/ε/P/E_elec
中间值和最终 ΔHf(模型)=−394.8679 kcal/mol；并说明真实分子从几何手算积分不可行，
可行路径是先取程序积分表（H、W、E_nuc）再手算 SCF。详见 §61。

---

## 25. 附录 A：可程序复算的计算条件（2026-10-25）

新增 H₂/PM7 完整算例：deck 条件、PM7 参数、积分表（S/H/γ/gss/E_nuc）、
化简 Fock 公式、两轮 SCF 迭代表、E_elec/E/ΔHf 计算及与 MOPAC 的逐项对照
（ΔHf=−31.74830 kcal/mol 完全一致）。见 `TRANSLATION_GAPS.md` §62。

---

## 26. 附录 B：积分公式（2026-10-25 增补）

新增 S 的闭式、H₁₂=βS、H₁₁=uss−V_AB、gss 参数、γ 的 reppd 公式、
E_nuc 的 ecore/reppd2 公式位置，以及程序验证步骤与精度目标；
H₂ 的 S/H₁₂ 可精确验证（<1e−10），γ/E_nuc 按完整实现可到 1e−6。
详见 `TRANSLATION_GAPS.md` §63。

---

## 27. 公式与验证汇编（2026-10-25）

新增 `FORMULAS_AND_VERIFICATION.docx`：单文档汇总全部公式、验证索引 V1–V17、
10 个手算案例、H₂ 计算条件与积分公式。详见 `TRANSLATION_GAPS.md` §64。

---

## 28. 附录 C 与完整参数表（2026-10-25）

新增 `MOPAC_parameters.xlsx` / `.csv`：MNDO、AM1、PM3、PM6、PM7、PM7-TS、RM1
七种方法、107 元素的全参数表；`FORMULAS_AND_VERIFICATION.docx` 增加附录 C
（方法对照、多方法 H₂ 参考值、任意元素验证步骤）。详见 §65。

---

## 29. 跨方法/分子验证矩阵（2026-10-25）

6 方法 × 7 分子 42 组 1SCF 对照：37 组支持组合的 HOF/E_elec/E_nuc 与官方
完全一致（含 FeO/PM6、PM7 的 d 轨道）；5 组因参数缺失与官方同样终止。
详见 `TRANSLATION_GAPS.md` §66 与文档附录 D。

---

## 30. 全部验证方案（2026-10-25）

新增附录 E：22 类分支状态矩阵、7 方法元素覆盖、10 项新增分支测试
（优化/COSMO/UHF/MOZYME/CI/ESP/PM7-TS/ZnO/LaO 全部 ΔHf=0）、
未完成项清单与“分支×方法×元素”全验证执行流程。详见 §67。

---

## 31. 积分公式全集与逐步验证（2026-10-25）

新增附录 F：重叠积分闭式/Gaussian 展开、单电子（H_μν=U+ΣV、βS）、双电子
（gss、γ=ev/√(r²+aee)、reppd 22 项、d 243 项）、核-核 ccrep 公式，
以及“几何→S→H→W→F→P→E→ΔHf”的 11 步验证清单与任意分子导出积分的方法。
详见 §68。

---

## 32. 验证顺序手册（2026-10-25）

新增附录 G：8 阶段验证顺序（符号→参数→积分→SCF→能量→方法×分子→分支→
记录/排查）、分支内部顺序、14 项勾选清单、失败排查顺序和扩展顺序；
首页增加①–⑧摘要。详见 §69。

---

## 33. 程序选项审计（2026-10-25）

按源码审计 305 个 keywrd 检查项，新增附录 H：分类、主要模式验证状态、
不可用/未验证/扩展项清单及“不能声称验证所有选项”的修订声明。
详见 §70。

---

## 34. 未完成分支审计（2026-10-25）

新增附录 I：stub/TODO 扫描 + 双程序实测，确认 FORCE/THERMO、FLEPO/DFP/
L-BFGS、STATIC、IRC、denrot、ENPART、MECI 分析、MOZYME+COSMO 等仍未完成/
异常；仅 EF 优化器通过。详见 §71。

---

## 35. 优化器修复（2026-10-25）

FLEPO/DFP 修复三个根因（hesinv 单位阵、ddot/supdot 1-based 偏移、linmin
指针偏移），水 PM7 结果与官方一致；L-BFGS 暂回退到 FLEPO。其余缺口见附录 I。
详见 §72。

---

## 36. FORCE 修复（2026-10-25）

fmat 改为 4 阶中心差分 Hessian；水 PM7 ZPE 10.118（官方 10.119），
频率基本一致。THERMO 仍崩溃、IR 强度未恢复。详见 §73。
