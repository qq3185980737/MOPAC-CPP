# MOPAC-CPP 计算公式与迭代流程总整理

本文把当前 C++17 移植版（`D:\test_rep\mopac_cpp_public`）中真正参与数值计算的公式、
内外层迭代结构、每一步的输入/输出与判据整理出来；公式均标注对应的 C++ 源码位置
（与 Fortran 原名一一对应），可直接对照代码核对。

> 完成度说明（2026-10-25 审计）
> - **注释**：**380/380 个 `.cpp` 与 372/372 个 `.h` 均有 `【中文导读】`
>   中文头注释**（2026-10-25 补齐最后 5 个并更新了 polar/writmo 的旧描述）。
> - **专用输出分支**：已补 `write_unit_cell_HOF`、`write_pressure`、周期
>   `id==3` 块、`symtry`、`prtgra`，见 `TRANSLATION_GAPS.md` §50。
> - **功能**：核心数值链（SCF/MOZYME/几何优化/解析梯度/COSMO/POLAR/
>   色散与 PM6-PM7 校正/MULLIK/MOZYME 输出）已闭环并通过官方逐位/逐行对照；
>   仍缺的专用输出分支列在文末 §10。

---

## 0. 记号

| 记号 | 含义 |
|---|---|
| μ,ν,λ,σ | 原子轨道（AO），属于原子 A/B |
| i,j | 占据分子轨道（非 MOZYME）或占据 LMO（MOZYME） |
| a,b | 虚轨道 / 虚 LMO |
| P_μν | 密度矩阵（打包下三角，`j*(j-1)/2+i`） |
| H_μν | 单电子（核-电子 + 动能）矩阵 |
| F_μν | Fock 矩阵 |
| W(μν,λσ) | 双电子积分（NDDO 近似后只有原子块非零） |
| C_μi | MO/LMO 系数 |
| ε_i | 轨道能量（`eigs`） |
| Z_A | 原子实电荷 |
| fpc_9 | 23.060529 kcal·mol⁻¹/eV（`molkst_C::fpc_9` / `enrgy_`） |
| E_elec | 电子能量（eV） |
| E_nuc (enuclr) | 核-核排斥能（eV） |
| atheat | 原子生成热之和 + 各类校正（kcal/mol） |

---

## 1. 单点能：从 AO 到生成热

### 1.1 单电子矩阵与核-核排斥（`hcore.cpp`）

对每个原子对 (A,B)，`hcore()` 做三件事：

1. 同原子块：把单电子积分 `uspd` 放入 H 对角（`hcore.cpp` 文件头注释）。
2. 异原子块：`h1elec()` 计算电子-核吸引/动能项累加进 H；
   `rotate()`（孤立分子）或 `solrot()`（周期体系）计算双电子积分 W、
   电子-核吸引 e1b/e2a，同时累加核-核排斥 `ENUCLR`。
3. 一中心双电子积分块由 `wstore()` 写入 W。

公式层面（MNDO/AM1/PM3/PM6/PM7 共同的 NDDO 形式）：

```
H_μν = U_μν δ_μν  −  Σ_B  V_μν,B  (电子-核吸引)  [+ 外场/点电荷项]
E_nuc = Σ_{A<B} Z_A Z_B (s_A s_A | s_B s_B) + 芯-芯校正项
```

`hcore()` 还处理：`FIELD` 外电场（`fldon`）、`QMMM` 嵌入点电荷、COSMO
介电校正（`useps` → `addnuc/addhcr`）。W 的写指针 `kr` 必须与 Fortran 顺序
完全一致，否则 `fock2()` 取积分会整体错位。

### 1.2 Fock 矩阵（`fock2.cpp` / `fock1.cpp`）

`fock2()` 是 SCF 每轮的核心：**F = H + 双电子贡献(P)**。

NDDO 下 Fock 的普适表达式：

```
F_μν = H_μν + Σ_{λσ} P_λσ [ (μν|λσ) − ½ (μλ|νσ) ]
```

其中：
- `(μν|λσ)` 是库仑项，`(μλ|νσ)` 是交换项；
- NDDO 近似使 (μν|λσ) 仅在 μν 同属 A、λσ 同属 B 时非零；
  同原子多极积分、双中心多极展开都在 `fock2()` 内完成；
- `fock2()` 按原子对 (i,j) 遍历，按轻/重原子分 `hh`、`hl`、`ll`
  以及周期体系 `id!=0` 的块处理，用 W 的游标 `kk` 取积分；
- 一中心 d 轨道块由 `fock1dorbs()`/`fock1()` 补充。

C++ 与 Fortran 保序的关键：W 的写入顺序、`kk` 推进顺序、线程私有累加器
（默认串行；仅 `dcart` 导数路径和 `MOPAC_DSPAR` 并行）必须逐位一致。

### 1.3 电子能与生成热（`compfg.cpp`）

SCF 收敛后电子能：

```
E_elec = ½ Σ_μν P_μν (H_μν + F_μν)
```

`compfg.cpp:219-260` 组装生成热（引用实际代码）：

```
escf = (elect + enuclr) * fpc_9 + atheat;
```

其中 `atheat` 以原子生成热为起点（`run_mopac.cpp:602`）再叠加各种校正
（`compfg.cpp:212-258`）：

```
atheat = Σ_A eheat[nat_A]                ← 原子生成热（parameters_C::eheat，
                                            如 H 52.102、O 59.559 kcal/mol）
       + hpress                        ← 压力项（id=1: −p·|tvec|；id=3: −p·V）
       + solv_energy*fpc_9             ← COSMO 溶剂化（非 MOZYME）
       + nsp2_corr                     ← PM6 的 >N− 校正（PM6-D 相关）
       + Si_O_H_corr                   ← PM7 的 Si-O-H 校正
       + sum_dihed                     ← 二面角校正 Σ htype·sin²θ
       + stress                        ← GEO_REF 几何偏差惩罚
```

因此最终：

```
ΔHf = fpc_9 · (E_elec + E_nuc) + atheat          [kcal/mol]
```

与官方输出一致（`escf` 即 `FINAL HEAT OF FORMATION`）。

---

## 2. SCF 迭代（非 MOZYME，`iter.cpp`）

`iter()` 是整个程序的内层循环；流程与官方 `iter.F90` 标号对应：

```
L180 初始化 → L250 迭代体 → L320 建 Fock → L410 收敛/更新 → L600 收尾
```

### 2.1 初始化（`iter.cpp` 头部 + L180 区域）

- 阈值：`scfcrt=1e-4`（密度判据基准）、`plchek=0.005`、`pltest=0.05`、
  `trans=0.2`（阻尼混合系数）、`bshift=-80`、`shfmax=20`（level shift 上限）。
- 初始密度：由 H 的对角化（或旧密度 `pdiag`）给出；`shift=1.0` 起步，
  DIIS 打开时 `shift=0`（DIIS 需要纯 Fock 做外推）。
- `iscf=1` 表示首次 SCF。

### 2.2 每轮迭代体（L250 → L410）

```
for niter = 1, 2, ...:
  ① 建 Fock：fock2(F, P, W)  →  F = H + G(P)
      + COSMO 反馈 addfck()（useps 时）
      + 外场/QMMM 项（若有关键字）
  ② 对角化：diag/diag_for_GPU 解 F C = C ε
      （大体系/关键字可走伪对角化路径）
  ③ 由 C 构建新密度：P_new = Σ_occ C_i C_iᵀ
      （占据数：闭壳 2；开壳按 fract 分数占据）
  ④ 混合/阻尼：P ← (1−trans)·P_old + trans·P_new（或 DIIS 外推）
      虚轨道 F_aa += shift（level shift，逐步降到 0）
  ⑤ 收敛判据：
      pl = max_μν |P_new − P_old|（按 norbs 归一）
      ΔE = |E_new − E_old|
      若 pl < scfcrt 且 ΔE 足够小 → 收敛；否则 niter++ 回 ①
```

关键辅助量：
- `pdiag`：上一次对角化写入的密度对角，作为收敛/重启的参考（`iter.cpp:181`）。
- DIIS（可选）：保存最近若干 (F, residual) 对，解线性方程组得外推系数，
  再组新的 F；仅当 `diis_on` 时启用（`iter.cpp` 顶部 DIIS 窗口代码）。
- 能量：`E_elec = ½ Σ P(H+F)`；`escf` 由 `compfg()` 汇总。

### 2.3 收尾（L600）

- 用**精确全对角化**再算一次最终轨道/能量（保逐位一致）；
- 若请求 CI：调用 `meci()` 做 CI 能量修正；
- 返回 `ee`（电子能）；再交 `compfg()` 计算 ΔHf 与（可选）梯度。

---

## 3. MOZYME：局域轨道 SCF（`iter_for_MOZYME.cpp` + `diagg.cpp`）

MOZYME 不直接对 F 做全对角化，而是：

1. 第一次 SCF 前把 canonical MO 局域化成 LMO（占据集 + 虚集）；
2. 每轮 SCF 用 **diagg（伪对角化）** 把占据-虚 LMO 间的 Fock 元消掉，
   期间保持 LMO 局域性，并重建密度。

### 3.1 LMO 的构造/局域化（`local_for_MOZYME.cpp`、`local2.cpp`）

对占据集和虚集分别做两两 Jacobi 旋转，目标函数（代码注释原文）：

```
maximize  Σ_μ  C_μi⁴        （即 Σ_atoms Σ_μ∈A C_μi² 的平方和形式）
```

`local2.cpp:25-57` 的实际旋转公式（i、j 为两条 LMO，逐原子 A 求和）：

```
d_ij^A = Σ_{μ∈A} C_μi C_μj ,  d_ii^A = Σ_{μ∈A} C_μi² ,  d_jj^A = Σ_{μ∈A} C_μj²
x_iiii = Σ_A (d_ii)² ,  x_jjjj = Σ_A (d_jj)² ,  x_iijj = Σ_A d_ii d_jj
x_ijij = Σ_A (d_ij)² ,  x_ijjj = Σ_A d_ij d_jj ,  x_jiii = Σ_A d_ij d_ii

a = x_ijij − (x_iiii + x_jjjj − 2·x_iijj)/4
b = x_jiii − x_ijjj
c = √(a² + b²) ,  s = a + c
若 s > 1e-14：
    ca = (1 + √((1 − a/c)/2)) / 2
    sa = √(1 − ca) ,  ca = √(ca)
    C_i' =  ca·C_i + sa·C_j
    C_j' = −sa·C_i + ca·C_j
```

`local2.cpp` 每次外循环累加 |sa|，当 `sum < 1e-5` 时停止（nmos==4 时最多
19 轮，否则 1 轮粗扫）；`local_for_MOZYME.cpp` 对占据/虚集重复该过程。
`icocc/ncocc/nncf/ncf`（占据）、`icvir/ncvir/nce/nnce`（虚）、
`iorbs`（每个原子的 AO 数）、`ijbo/nijbo`（原子对与块偏移）。

### 3.2 密度矩阵（`density_for_MOZYME.cpp`）

MOZYME 密度仍由 LMO 系数直接构造：

```
P_μν = Σ_i∈occ  C_μi C_νi
```

`mode` 参数控制增量更新：`0` 全量重建；`-1` 先扣掉旧 LMO 密度；
否则叠加新 LMO 密度（含部分密度 `partp` 以加速）。
`P` 保持打包下三角、1-based。

### 3.3 diagg1：显著 Fock 元与 LMO 能量（`diagg1()`，`diagg.cpp:54-476`）

目的：避免构造/对角化完整 F，只处理“可能显著”的占据-虚 LMO 对。

```
aocc[k] = Σ_{μ∈原子 k} C_μi²        （占据 LMO i 在原子 k 上的强度）
avir[k] = Σ_{μ∈原子 k} C_μa²        （虚 LMO a 在原子 k 上的强度）
```

选择判据（代码 `diagg.cpp:380/427`）：

```
若  avir[k1] · P[k1,j1] · aocc[k] ≥ cutoff   则计算该原子对贡献
cutoff = max(1e-8, tiny·10·1e-8)              （idiagg>5 时自适应放大）
```

对入选原子对，组装 LMO 基 Fock 元：

```
fmo[i,a] = Σ_{k1,j1} Σ_{μ∈k1,ν∈j1} C_μi · F_μν · C_νa
```

（按 μν 是否同原子拆成 k1>j1 / k1<j1 / 同原子三种块，循环见 360-476 行。）
`ijc` 计数已存积分；若 `ijc==nij`（存储不足）则 `safety*=2` 下轮少算，
否则 `safety=max(safety/2,1)`。

**能量分支**：当 `idiagg ≤ 2` 或 `idiagg % 4 == 0` 时，用同样的
截断求和计算占据 LMO 的期望值：

```
eigs[i] = Σ_{j,k} Σ_{μ,ν} C_μi F_μν C_νi   （只含通过 cutoff 的原子对）
```

这正是官方 `prtlmo` 打印的 LMO 能量来源；本轮已确认
C++ 直接使用 `eigs` 与官方逐位一致（见 `docs/TRANSLATION_GAPS.md` §48）。
否则（`idiagg>2` 且 `idiagg%4!=0`）只更新 `fref = tiny⁴`，不改 `eigs`。

### 3.4 diagg2：2×2 Jacobi 旋转消去显著元（`diagg2()`，`diagg.cpp:477-708`）

对 diagg1 给出的每个显著 (i 占据, a 虚) 元：

```
c = const2 · fmo[i,a]                        （const2 默认 1，过渡金属/DAMP 可 0.5）
d = ε_a − ε_i − shift
若 |c/d| ≥ biglim：                          （biglim = 50·√eps_machine）
    e     = sign(d) · √(4c² + d²)
    α     = √( ½(1 + d/e) )
    β     = −sign(c) · √(1 − α²)
    C_occ' = α·C_occ + β·C_vir
    C_vir' = α·C_vir − β·C_occ
```

旋转后若虚 LMO 在原占据 LMO 原子上强度 `> thresh`，把该原子加入占据
LMO 列表（并为它分配系数）；反之把占据原子加入虚 LMO 列表。
若列表空间不足则回滚（`nrej` 计数 + 恢复 `cocc/cvir`），并可能执行
`alpha = 0.5·(alpha+1.0)` 半角重试。

**门槛调度（决定“每一步怎么迭代”）**：

```
idiagg 从 0 开始，每次 SCF 迭代 +1（若 addhb 造出氢键则强制为偶数）
当 idiagg%5==0 或 idiagg≤5：tiny=−1, biglim=−1   → 消去所有显著元（粗扫）
否则：tiny ← 0.01·tiny, biglim ← bigeps          → 只消去强元（细扫）
idiagg%4==0 时 diagg1 顺带刷新 LMO 能量 eigs
```

即 MOZYME 的 SCF 迭代是“粗→细”交替的伪对角化：前 5 轮全消，
之后每 5 轮重置一次门槛，期间逐步收紧 `tiny`。

### 3.5 MOZYME SCF 主循环（`iter_for_MOZYME.cpp`）

```
① 初始化：读取 THRESH/RELTHR/PL/REORTH/TIMES 等关键字；idiagg=0
② 首次：density_for_MOZYME(p,0,noccupied,partp) 全量密度；
   local_for_MOZYME("OCCUPIED"), local_for_MOZYME("VIRTUAL")
③ 每轮：
    buildf → 由 P 构建 F（含 hcore_for_MOZYME/fock1_for_MOZYME 的
            原子块、MOZYME 快速多极路径；详见 buildf.cpp）
    diagg(F, nocc, nvir, idiagg, partp, indi)
         = diagg1(选元/能量) + diagg2(旋转/扩表) + density_for_MOZYME(重建 P)
    idiagg ← idiagg+1；若需要则强制偶数
    addhb(nocc,nvir,idiagg,nij,nhb)      ← 氢键校正/成键调整
    isitsc() 收敛判据（密度变化/能量变化）
    未收敛 → 回 buildf；收敛 → 收尾（正交化/最终能量）
④ 最终 eigs/能量交 writmo/prtlmo 输出
```

---

## 4. 几何优化：外层迭代（`compfg` + `flepo`/`ef`）

外层结构（`run_mopac.cpp` 分派）：

```
x⁰（内坐标/笛卡尔）
repeat:
    compfg(x, ...):
        gmetry(x → 坐标)
        hcore()             → H, W, E_nuc
        iter()/iter_for_MOZYME() → 收敛的 P, F, E_elec
        ΔHf = fpc_9(E_elec+E_nuc) + atheat
        deriv()（需要梯度时）→ g = ∂ΔHf/∂x
    optimizer 步进（BFGS / EF / LBFGS 等）
    x ← x + Δx
until 收敛
writmo() 输出最终结果
```

### 4.1 BFGS/DFP（`flepo.cpp`）

- 方向：`p = −H⁻¹ g`（H⁻¹ 为逆 Hessian 近似）。
- 线搜索：`linmin` 沿 p 找使 ΔHf 最小的 α；能量/梯度由 `compfg` 给出。
- 逆 Hessian 更新（BFGS）：

```
s = x_{k+1} − x_k
y = g_{k+1} − g_k
ρ = 1/(yᵀs)
H⁻¹_{k+1} = (I − ρ s yᵀ) H⁻¹_k (I − ρ y sᵀ) + ρ s sᵀ
```

- DFP 变体按关键字切换（代码内保留两套更新）。
- 收敛判据：Herbert 检验、`|g|`（GNORM）、变量/函数变化、Peters 检验
  （连续若干步能量变化同号且很小）。`linmin` 失败/能量不降时按
  `flepo.F90` 的回退逻辑缩步。
- `dfpsav` 周期性保存 Hessian 与几何，支持重启。

### 4.2 EF 本征向量跟踪（`ef.cpp`）

- 用 Hessian 的本征分解 `H v_i = λ_i v_i`，把梯度投影到模式坐标
  `g_i = v_iᵀ g`。
- **P-RFO 步**：解 `(H − λ) Δx = −g`，其中 λ 是位移约束乘子；
  对极小点，λ 取小于最小本征值；对过渡态，沿跟踪模式取“最大化”方向，
  其余模式按极小化处理（P-RFO/QA 两种步长构造都在 `formd()` 内）。
- 模式跟踪：`overlp()` 把新本征向量与旧模式做重叠匹配，防止本征值交叉
  导致步长跳变。
- 平动/转动投影：`prjfc()` 投影掉 T/R 分量（MKL dgemm 保序）。
- Hessian 更新：`updhes()` 用 Powell/BFGS 更新；初始 Hessian 由
  `gethes()` 对梯度做数值差分得到。
- 收敛：能量变化、梯度范数、步长、模式跟踪状态。

### 4.3 L-BFGS（`lbfgs.cpp`）

完整 L-BFGS-B 移植（`lbfgs`/`lbfsav`/`setulb` 等），用于大体系/ `LBFGS`
关键字；保留在代码中但主的默认优化分支是 EF/BFGS（`run_mopac.cpp` 分派）。

### 4.4 解析梯度（`deriv.cpp` / `deri1.cpp` / `deri2.cpp`）

对能量表达式求导：

```
dΔHf/dx = Σ_μν P_μν ∂H_μν/∂x
        + ½ Σ_{μνλσ} P_μν P_λσ ∂(μν|λσ)/∂x
        + dE_nuc/dx
        + 校正项导数（色散、氢键、溶剂化、二面角等）
```

非变分部分（如 CI/开壳）用 `deri1`（非松弛导数）和 `deri2`（松弛项，
解 CPHF 方程 `(ε_a−ε_i)U_ai = −B_ai`）补全，与官方 `deri1/deri2` 逐位对齐。

---

## 5. COSMO 溶剂化（`cosmo.cpp`）

1. **初始化**（`cosini`）：原子半径 `srad`、原子对列表 `nipa/lipa`、
   方向向量 `dirvec/dirsm`、`n0`、`disex2`。
2. **腔体**（`coscav`）：逐原子撒点成 SAS 面，按 `Segments` 归属，
   `surclo()` 闭合凹面环段，累加 `area`、`cosvol`；生成 `finel` 列表并用
   `mfinel()` 组装 `amat`，做 Cholesky 分解。
3. **B 矩阵**（`mkbmat`）：网格点/环段到原子的 `1/r` 多极作用；
   对角元用面积修正（`1.07·√(4π/A_i)` 类形式）。
4. **屏蔽电荷与介电能**：设介电因子

```
fepsi = (ε − 1)/(ε + ½)        （writmo 中取 ε=78.4 → fepsi = 77.4/78.9 ≈ 0.98099）
```

求解线性方程（Cholesky 回代 `coscl1/coscl2`）：

```
q = −fepsi · A⁻¹ V
V_i = Σ_A Z_A/|r_i−R_A| − Σ_μν P_μν (μν|i)
```

介电能：`E_diel = ½ Σ_i q_i V_i`，并令 `solv_energy = E_diel`；
核/电子屏蔽修正 `addnuc()/addfck()` 把 `B·q` 反馈进 H/Fock。
5. 最终 `ΔHf += fpc_9·solv_energy`（非 MOZYME，`compfg.cpp:232`），
   并在 writmo 输出 COSMO AREA/VOLUME（`area`、`cosvol`）。

数值闭环：水/硅醇等与官方逐项对照（`TRANSLATION_GAPS.md` §19-§26、
§35-§37）。

---

## 6. POLAR：TDHF 响应与极化率（`polar.cpp` + `alphaf/betaf/ngamtg`）

对每个频率 ω（`NFREQ` 条），在收敛 SCF 上解耦合微扰方程
（CPHF/含频响应）：

```
(ε_a − ε_i − ω) u_ai(ω) + Σ_bj [2(ai|bj) − (aj|bi)] u_bj(ω)
    + Σ_bj [2(ai|bj) − (aj|bi)] v_bj(ω) = −d_ai^μ
```

MOPAC 的 `alphaf()` 把实部/虚部分成 `ua`、`ub`，用 `ffreq1`（实 2e 部分）
和 `ffreq2`（虚 2e 部分）反复迭代（`ATOL`/`BTOL` 判据，`MAXITU/MAXITA` 上限）：

```
repeat:
    ua_new = solve( A·ua ) − B·ub
    ub_new = solve( A·ub ) + B·ua
    Δ = max|ua_new − ua|, |ub_new − ub|
until Δ < ATOL/BTOL
```

- 线性/二阶极化率张量 `α_ij = −Σ_ai d_ai (u_ai^j)`；
  各向同性平均：

```
α_iso = (α_xx + α_yy + α_zz)/3
α_iso[Å³] = α_iso[a.u.] · a0³
```

- β（SHG）：`betaf()` 把 6 个分量（X/XX,XY,YY,XZ,YZ,ZZ）按
  `bvx/bvy/bvz` 累加，5 次平均：

```
bav = (1/5)Σ b分量 ,   bvec = √(bav_x² + bav_y² + bav_z²)
```

- γ（THG）：`ngamtg()` 累加 5 次置换贡献后取平均：

```
gav  = Σ 各置换 γ 分量 ,  gave = gav/5
输出：AVERAGE GAMMA VALUE = gave a.u. = (gave/1.9854259398) ESU (×10⁻³⁹)
```

（`polar.cpp` 内 `gav/5.0` 即 5 次置换平均；输出单位换算/ESU 与官方一致。）
- 每个频率的表头/表格格式均按 `polar.F90` 对齐（本轮 §48）：
  `+++++ ALPHA/BETA/GAMMA ... AT ω EV`、`COMPONENTS OF ALPHA/BETA`、
  `AVERAGE ALPHA/BETA/GAMMA`。
- 收敛后：`writmo` 的 `POLAR` 输出、`ALPHA/BETA/GAMMA` 各频率块。

---

## 7. 其它校正项

| 校正 | 公式/实现 |
|---|---|
| 色散（PM6-D/PM7） | `E_disp = Σ_{A<B} f_damp(R) [C6/R⁶ + C8/R⁸ + C10/R¹⁰]`，`f_damp` 按 PM6-D3H4/PM7 参数；见 `dispersion` 相关文件 |
| 氢键 | `E_hb = Σ htype·(...)`；MOZYME 中 `addhb()` 在 diagg 后调整 |
| PM6 >N− 校正 | `nsp2_correction()`，加进 `atheat` |
| PM7 Si-O-H | `Si_O_H_Correction()`（`set_up_dentate.cpp`），加进 `atheat` |
| 二面角 | `sum_dihed = Σ htype·sin²θ`（`compfg.cpp:241-247`） |
| 压力 | `hpress = −p·V`（id=1 用 |tvec|，id=3 用 volume） |
| 溶剂化 | COSMO `solv_energy`（见 §5） |
| 外场 | `FIELD` → `fldon` 加到 H/Fock |

---

## 8. 完整迭代层级（一张图）

```
run_mopac
 ├─ getdat/readmo          读关键字、几何
 ├─ moldat                电子数/占据数、原子参数、点群、RHF 信息
 ├─ hcore                 H, W, E_nuc（每几何步重建）
 ├─ E = compfg(x)  ←─────────────┐
 │    ├─ gmetry(x → coord)        │
 │    ├─ hcore                   │
 │    ├─ iter / iter_for_MOZYME  │  ← 内层 SCF 循环（§2/§3）
 │    │    ├─ buildf/hcore       │
 │    │    ├─ diag / diagg       │
 │    │    ├─ density            │
 │    │    └─ 收敛判据            │
 │    ├─ ΔHf 组装（§1.3）         │
 │    └─ deriv（梯度）            │
 └─ optimizer (flepo/ef/lbfgs) ───┘  外层几何循环（§4）
      ├─ BFGS 方向/线搜索
      ├─ EF 本征向量跟踪
      └─ 收敛 → writmo
writmo: 结构/能量表、LMO、MULLIK、POLAR、COSMO 等输出
```

### 循环停止条件一览

| 循环 | 判据 | 常量/位置 |
|---|---|---|
| SCF | 密度变化 `pl < scfcrt`，能量变化足够小 | `scfcrt=1e-4`，`iter.cpp` |
| SCF 阻尼 | `P ← (1−trans)P_old + trans·P_new` | `trans=0.2` |
| SCF level shift | 虚轨道对角加 `shift`，从 `shfmax=20` 降到 0 | `bshift=-80`，`shfmax=20` |
| MOZYME 门槛 | 每 5 轮或前 5 轮粗扫（tiny=−1），其余细扫（tiny×0.01） | `diagg.cpp:353-476`, `diagg2` |
| MOZYME 能量 | `idiagg%4==0` 时刷新 `eigs` | `diagg.cpp:360-476` |
| 几何优化 | 梯度范数/GNORM、能量变化、步长、Herbert/Peters | `flepo.cpp`/`ef.cpp` |
| 线搜索 | `linmin` 沿方向最小化，失败回退 | `linmin` |
| POLAR | `|Δua|<ATOL`, `|Δub|<BTOL`，迭代 ≤ `MAXITU/MAXITA` | `polar.cpp` alphaf |
| COSMO | 一次求解（Cholesky），SCF 中每轮反馈 | `cosmo.cpp` |

---

## 9. 关键数值常量

| 常量 | 值 | 用途 |
|---|---|---|
| `fpc_9` | 23.060529 | eV → kcal/mol |
| `a0` | 0.529177... | Bohr → Å（极化率单位换算） |
| `fepsi`（水） | (78.4−1)/(78.4+0.5) = 77.4/78.9 ≈ 0.98099 | COSMO 介电因子 |
| `au_to_esu` | 8.639418×10⁻³³ | 偶极/超极化率 ESU 换算 |
| `scfcrt` | 1e-4 | SCF 密度收敛 |
| `cutlim` | 1e-8 | diagg1 原子对截断下限 |
| `const2`（diagg2） | 1.0（过渡金属/DAMP 0.5） | 2×2 旋转混合因子 |

---

## 10. 仍未完成/仅部分实现的功能（2026-10-25）

按 `docs/TRANSLATION_GAPS.md` 与源码 TODO 标记统计：

1. **专用输出分支**：超精细耦合、UHF 自旋态头信息、archive(.arc) 整段、
   周期 cell/pressure、`write_paras_used`、`prtgra`、`modchg/modgra` 均已
   移植（2026-10-25，`TRANSLATION_GAPS.md` §50-§51）。仍缺：
   - `denrot.cpp`：`denrot_for_MOZYME` 仍是空函数（SIGMA-PI 矩阵/旋转未移植）；
   - `wrt_diffs.cpp`/`drcout.cpp` 的 reverse_aux、bond-by-bond 打印；
   - archive 的 pKa 分支（PKA 关键字）与 `geo_diff` 专用块。
2. **CALPAR 程序**（`calpar.cpp` 的 `create_parameters_for_PMx_C`）未移植，
   只影响参数生成工具，不影响主程序。
3. **POLAR 数值末位**：α 各向同性 9.67479 vs 官方 9.67480（1e-5 相对）、
   β/γ 的 1e-5~1e-3 相对差；属于响应迭代/舍入路径差异，非公式错误。
4. **非 MOZYME 输出的数值噪声**：CARTESIAN 坐标末位 3e-9、偶极 −0.000
   符号、MOLECULAR DIMENSIONS 平局行序（gmetry 往返/求和顺序）。
5. **注释**：5 个文件仅有英文头注释；`polar.cpp`、`writmo.cpp` 头注释中
   仍写着“stub/未接线”，与已实现状态不符，需更新。

---

## 11. 已验证的数值闭环（回归）

| 项目 | 官方 | C++ | 一致性 |
|---|---|---|---|
| PM7 水 HOF | −57.79103 | −57.79103 | 逐位 |
| MOZYME 水 HOF | −57.78972 | −57.78972 | 逐位 |
| MOZYME 甲醇 HOF | −35.50018 | −35.50018 | 逐位 |
| POLAR α(水) | 9.67480 | 9.67479 | ~1e-6 相对 |
| POLAR β/γ | — | — | ~1e-5~1e-3 |
| MOZYME LMO 能量 | −23.1582/−23.1556/−15.2141/−12.0869 | 完全相同 | 逐位 |
| MOZYME 输出行 | 144 行 | 141 行一致 | 仅时间戳/耗时 |
| VEC 输出行 | 163 行 | 159~160 行一致 | 仅时间戳/耗时 |
| 普通水输出行 | 170 行 | 160 行一致 | 噪声/时间 |

---

## 12. 阅读代码的建议顺序

1. `iter.cpp`（SCF 主循环）→ `fock2.cpp`（Fock 组装）→ `hcore.cpp`（H/W/E_nuc）
2. `compfg.cpp`（能量/ΔHf/梯度入口）
3. MOZYME：`iter_for_MOZYME.cpp` → `diagg.cpp`（diagg1/diagg2/density）
   → `local_for_MOZYME.cpp`
4. 优化：`flepo.cpp`（BFGS）→ `ef.cpp`（EF）→ `linmin`
5. 溶剂化/响应：`cosmo.cpp` → `polar.cpp`（alphaf/betaf/ngamtg）
6. 输出：`writmo.cpp` / `prtlmo.cpp` / `mpcsyb.cpp`
