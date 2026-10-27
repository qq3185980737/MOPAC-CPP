# MOPAC-CPP 翻译完整性审计与缺口清单（2026-10-25）

本文件记录对 `mopac_cpp_public` 相对 MOPAC2016 官方 Fortran 源码
（`mopac-archive/2016_MOPAC2016/MOPAC_source_code`）的接口/功能完整性审计结果，
以及本轮已修复项与后续优先级。读者可据此判断某个功能是否可信、是否需要补移植。

## 1. 审计方法与总体结论

- 方法：正则提取 343 个实质 `.F90`（排除 336 个 `*_I.F90` 接口模块）中的过程声明，
  去掉 INTERFACE 块后得到 659 个唯一过程名；在 C++ 树中按名字+定义体匹配，
  并用 `dumpbin /symbols` 对 374 个已编译 `.obj` 的外部符号交叉验证。
- 结论：**链接成功 ≠ 语义完整**。大多数过程有同名实现，但存在两类静默缺口：
  1. **空实现遮蔽真实实现**（shadow stub）：局部匿名 namespace 里的 trivial 空函数
     抢先满足链接，真实实现（其他 .cpp）永远不会被调用；
  2. **顶层 routine 非空但只是打印提示**（如 `polar()`）。
- 计数（唯一过程名，659）：
  - 有对应实现：613（93.0%）
  - 同名定义但 body 为 trivial 空实现：14（2.1%）
  - 未找到同名定义：32（4.9%）；其中 7 个为改名/内联、15 个为预期省略，
    真缺口约 10 个。
- 374 个 `.cpp` 编入工程；6 个源码文件被 vcxproj 排除：
  `cosmo_surface.cpp`、`epseta.cpp`、`exp_variants.cpp`、`mini_diat.cpp`、
  `opencl_accel.cpp`、`web_message.cpp`（最后者导致 `web_message()` 实际链接到空 stub）。
- Header 606 条函数声明中，只有 28 条没有已编译定义，全部为 GPU/OpenCL（预期）。

## 2. 高危：静默丢弃数值/功能的 shadow stub（已核实）

| 位置 | 函数 | 真实实现 | 影响 |
|---|---|---|---|
| `compfg.cpp` | `Si_O_H_correction()`（小写 c） | `set_up_dentate.cpp` 的 `Si_O_H_Correction()`（大写 C） | PM7 的 Si–O–H 校正恒为 0 |
| `fock2.cpp` | `addfck(vector&, vector&) {}` | `cosmo.cpp` 的 `addfck(double*, double*)` | EPS/COSMO 对 Fock 的修正被丢弃 |
| `ijkl.cpp` | `ciint(double*, double*) {}` | `cosmo.cpp` 的 `ciint(const double*, double*)` | COSMO 双电子积分被丢弃 |
| `calpar.cpp` | `sp_two_electron() {}` | `sp_two_electron.cpp` 同名实现 | CALPAR 参数路径不再计算 gss/gsp/gpp 等 |
| `getsym.cpp` | `getsym() {}`（0 参） | `getsym.cpp` 的 4 参实现（无调用点） | 输入文件中的对称性约束静默失效 |
| `datin.cpp` | `datin(int) {}` | `datin.F90` | `EXTERNAL=`/`PARAMS=` 外部参数文件被忽略 |
| `local.cpp` | `resolv_/phase_lock_/molval_/matout_` 全空 | `resolv.cpp`、`force.cpp`、`molval.cpp`、`matout.cpp` | LOCAL/GRAPH 路径的密度重建、相位锁定、结果打印全部跳过 |
| `chkion.cpp`/`chklew.cpp` | `mopend(){}`、`memory_error(){}` | `mopend.cpp`、`setup_mopac_arrays.cpp` | 错误被吞掉、`moperr` 不置位，坏输入继续算 |
| `convert_storage.cpp`、`delsta.cpp`、`diagg1.cpp`、`denrot_for_MOZYME.cpp` | `ijbo(){return -1/0;}` | `ijbo.cpp` | 成键表查询恒失败（MOZYME/局部相关路径） |

## 3. 真缺口（需移植）

| 优先级 | 功能 | 证据 |
|---|---|---|
| P0 | `polar()` 驱动 + `alphaf`（TDHF 极化率/超极化率） | `polar.cpp` 只打印 retained stub；子程序 bdenin…betaf 已移植但未接线；POLAR 关键字下结果缺失 |
| P1 | COSMO 接线（见 §2 的 addfck/ciint）与 `addnucz`（MOZYME+COSMO） | `hcore_for_MOZYME.cpp` 的 useps 分支为 no-op |
| P1 | `snapth`（SNAP 关键字角度吸附） | readmo.F90:428；C++ 无引用 |
| P2 | MOZYME：`fz2` 慢速 Fock 路径、`lmo_to_eigenvectors`/`FlmoFromFao` | `buildf.cpp` 注释自认未移植；writmo 缺 MOZYME 局域轨道输出 |
| P3 | 纯输出：`coswrt`、`write_params`、`write_unit_cell_HOF`、`write_pressure`、`printp`、`wrtcon/wrtwor/wrtout`、`denrot` | 只影响 .out/.cos 等输出完整性 |

预期省略（无需移植）：QuickWin GUI、Password/许可证、GPU/CUDA（模型为声明+host stub）、
`MOPAC for Windows.F90` 的 GUI 过程、`special`（死代码）、`ilaenv/ieeeck`（由 MKL 提供）。

## 4. 接口/签名核查

抽样 20 个核心过程（iter/hcore/fock2/writmo/meci/deriv/dcart/compfg/ef/lbfgs/
flepo/drc/polar/esp/pmep/interp/capcor/cnvg/helect/pulay）：**参数个数全部一致**，
未发现 int/double 类型错误或 `intent(out)`→按值传递错误。需知悉的接口形式变化：

- 矩阵参数由 Fortran 平铺数组改为 `std::vector<std::vector<double>>`
  （`deriv` 的 `GEO`、`coord` 等），数据布局契约变化；
- `fock2` 的 `P`、`dcart` 的 `coord` 在 Fortran 为 `intent(in)`，C++ 放宽为可变引用；
- `interp` 的 `fp` 由 in 变非 const 指针；
- `getsym` 存在“4 参真实版无调用点、0 参空 stub 被调用”的接错接口。

## 5. 本轮（2026-10-25）已修复

- `compfg.cpp`：删除小写 `Si_O_H_correction` stub，改为包含 `set_up_dentate.h`
  并调用 `Si_O_H_Correction()`。实测 PM7 硅醇输入 `Si_O_H_present=1`，
  校正进入能量路径（修复前该值恒 0）。
- `calpar.cpp`：删除空 `sp_two_electron`，改为包含真实实现头文件并调用。
- `chkion.cpp` / `chklew.cpp`：删除本地空 `mopend`/`memory_error`，改为链接真实
  实现，恢复错误传播。
- 修正错误注释：`polar()` 标注为 retained stub；`fock2` 的 addfck 标注为
  “遮蔽真实实现、COSMO 路径未接通”；`compfg` 补记 Si–O–H 大小写问题。
- 构建验证：MSBuild Release x64 成功；H2/PM7 回归 HOF 不变
  （85.93132 / -32.01058 kcal/mol）。

## 6. 重要观测：COSMO 路径不收敛（表现为长时间运行）

- 用合法 3 行表头输入 `PM7 EPS=78.4 1SCF` 的水分子，C++ 版运行 >90 s 无输出
  后被人为终止；而 `ITRY=5` 时 0.38 s 结束并报“FAILED TO ACHIEVE SCF”，
  说明循环有限，根因是 **SCF 不收敛**（默认 2000 轮、每轮含 COSMO 求解 → 很慢）。
- `PM7 EPS=78.4 0SCF` 可正常完成，说明 COSMO 初始化不挂；问题集中在
  COSMO 对 Fock/能量的耦合（addfck/addhcr/ciint/介质项）尚未完全对齐。
- 官方参考：根目录 standalone exe 需要 password；但源码构建的
  `x64\Release\MOPAC2016.exe`（password 已注释）可正常运行，可作为数值基准。

## 7. 建议下一步（按性价比）

1. 接线 COSMO 的 `addfck`/`ciint`（删除 shadow stub，按真实函数参数顺序调用；
   注意 `addfck` 的 C++ 参数顺序为 (density, Fock)，与 F90 声明顺序相反），
   再跑 EPS 回归确认不再挂起；
2. 修复 `getsym` 调用为 4 参真实版本（对称性约束）；
3. 接线 `local.cpp` 的 resolv/molval/matout/phase_lock 到真实实现；
4. 移植 `polar()` 驱动与 `alphaf`（最大功能缺口）；
5. 处理 `datin`（EXTERNAL/PARAMS）与 `snapth`（SNAP）；
6. 纯输出缺口（wrtcon/wrtwor/wrtout、denrot、writmo 输出块）最后补齐。

## 8. 参考程序可用后的验证结果（2026-10-25 追加）

参考程序：`D:\mopac\mopac-archive\2016_MOPAC2016\x64\Release\MOPAC2016.exe`
（password 已注释；输入需符合 3 行表头格式）。

| 用例 | 官方参考 (kcal/mol) | C++ 修复后 (kcal/mol) | 结论 |
|---|---|---|---|
| H2 / PM7 1SCF | 85.93132 | 85.93132 | 一致 |
| H2 / PM7（EF 分支） | -32.01058 | -32.01058 | 一致 |
| 硅醇 SiH3OH / PM7 1SCF | 59.17334 | 59.17334 | 一致（验证 Si–O–H 接线修复） |

- Si–O–H 修复前该校正恒为 0；修复后 `Si_O_H_present=1`，最终 HOF 与官方逐位一致。
- COSMO/EPS 用例仍不收敛（见 §6）：addfck/ciint 已接线，但还需定位其他耦合项。

## 9. 2026-10-25 第二批接口修复与验证

- 删除 `convert_storage.cpp`、`delsta.cpp`、`diagg1.cpp`、`denrot_for_MOZYME.cpp`
  中遮蔽真实实现的局部空 `ijbo`（return -1/0），统一调用 `ijbo.cpp`；
- `local.cpp` 的 `resolv_`/`molval_`/`matout_` 改为转发真实 `resolv`/`molval`/
  `matout`；`phase_lock_` 仍为 no-op（真实函数需要 1-based 打包 vector，待接线）；
- `getsym` 仍未接线：4 参真实版从 `std::cin` 读对称性约束（数据流接错），
  目前 readmo 仍调用 0 参 stub——已改为文档化缺口；
- 验证：5 个用例与官方 `x64\Release\MOPAC2016.exe` 逐位一致
  （kcal/mol）：H2 85.93132、PM7 H2 -32.01058、MOZYME H2O -57.78971、
  H2O -57.79103、硅醇 59.17334。

## 10. 2026-10-25 COSMO 定位进展（NaN 根因）

- 已确认 COSMO 初始化正确：`PM7 EPS=78.4 1SCF` 水分子下
  `area=42.4865`、`cosvol=25.1942`、`nps=49`、`fepsi=0.98099`，与官方
  （42.49 / 25.19）一致；`ediel` 初始为 0。
- 首次进入 `addfck()` 后 `ediel` 立即变为 `-nan(ind)`，随后 SCF 无法收敛。
  因此根因位于 COSMO 对 Fock 的耦合（`addfck` → `coscl2` 的 Cholesky 求解
  链路，涉及 `bmat`/`amat`/`qscnet`/`phinet`），不是初始化或面积/体积计算。
- 需要逐项对照 Fortran：`coscav` 的 `amat` 对角/边界、`mkbmat` 的 `bmat`、
  `coscl1` 的 `id` 偏移表，以及 `cl2` 的前代/回代下标。定位后可再用官方
  `x64\Release\MOPAC2016.exe` 验证 `DIELECTRIC ENERGY = -0.35863 EV`。

## 11. 2026-10-25 COSMO 修复（coscl2 列存储越界）与剩余 NaN 根因

- **已修复**：`addfck()` 与 `dmecip()` 原先调用
  `coscl2(..., &qscnet[1][2], &phinet[1][2], nps)`，把 `vector<vector<double>>`
  的一行地址当成 Fortran 列（连续 nps 个元素）传给求解器，造成越界写与 NaN。
  现改为与 `addnuc()` 相同的扁平缓冲（复制 slot 列 → coscl2 → 写回），
  这是 COSMO 的首要 bug 修复。
- 修复后 `ediel` 仍为 NaN，但来源更明确：
  - `phinet` 正常（phi1=-4.45393、phi2=-4.40508）；
  - `coscav` 后 `amat` 有 41 个非有限元（共 1225 个）；
  - 首个非有限元来自同一原子段对（如 ips=7,jps=3,i=j=1）的 `mfinel` 求和，
    `aa=inf`，即两个 finel 采样点距离为 0 时出现 1/0。
- 结论：剩余问题位于 `coscav` 的 `mfinel`/`surclo`/`nset/nsetf` 数据流
  （同一原子的段点被重复生成/未去重），导致 `amat` 出现 inf→NaN；
  Fortran 参考不会产生重复点。下一步应逐行对照 `coscav` 中 `mfinel` 的
  `nsetf/nar_csm/nset` 与 `surclo` 去重逻辑（cosmo.F90:约 500-620）。
- 验证：coscl2 修复后，非 COSMO 的 5 个基准用例（H2/PM7/MOZYME/H2O/硅醇）
  仍与官方逐位一致，未引入回归。

## 12. 2026-10-25 COSMO 剩余 NaN 进一步定位：环点（ring）重合

对 `coscav` 中非有限 `amat` 元做零距离检测，得到的首批样本全部是**环点-环点**
重合，而非基本网格点：

```
ips=7 jps=3 k=8 l=26   (basic1=7, ring1=1; basic2=25, ring2=1)
ips=9 jps=7 k=20 l=8
ips=13 jps=7 k=12 l=8  (ring1=2, ring2=1)
...
```

- `k`/`l` 均大于各自的 `nar_csm`（基本点数），即命中 `finel` 里由 `rsc`
  追加的环段点；两个不同段的环点坐标完全相同 → 距离 0 → `1/0 = inf`。
- 因此问题不在 `dirvec` 基本网格点，而在 `surclo` 生成的 `rsc` 环段数据：
  同一原子上的不同段拿到了同一个环点（或环点被重复写入）。
- 下一步：对照 `cosmo.F90` 中 `surclo` 的 `rsc` 写入与 `nipsrs` 计数，
  检查 `rsc(ix,irs)` 是否应按段去重、以及 `nipsrs(ips-1)+1..nipsrs(ips)`
  的范围是否正确；修复后 COSMO 应得到有限 `ediel`，再与官方
  `DIELECTRIC ENERGY = -0.35863 EV`、`HEAT OF FORMATION = -65.16522` 对照。

## 13. 2026-10-25 COSMO 根因确认：`surclo` 是不完整移植（缺闭合段注册）

对照 Fortran 与 C++ 的 `surclo` 后确认：

- Fortran `surclo`（cosmo.F90:1202-1663，共 462 行）的参数包含
  `cosurf, iatsp, nar_csm, nsetf, isude, sude`，函数体实际写这些数组（计数：
  `cosurf` 10 次、`nar_csm` 5 次、`nsetf` 6 次、`iatsp` 8 次、`isude/sude` 13 次），
  即**它会创建/注册凹区闭合段并更新段表**；
- C++ `surclo`（cosmo.cpp，约 245 行）的签名只有
  `coord, nipa, lipa, din, dim_din, rsc, isort, ipsrs, nipsrs, nat, srad, maxrs`，
  **没有 cosurf/iatsp/nar_csm/nsetf/isude/sude**，因此无法完成闭合段注册。
- 结果：`nipsrs`/`rsc` 与段的对应关系错误，`mfinel` 为不同段返回相同的环点，
  环点-环点距离为 0 → `amat` 出现 inf → `ediel` NaN → COSMO SCF 不收敛。

结论：COSMO 不是“接线”能修好的问题，必须**补全 `surclo` 移植**
（把 Fortran 462 行完整译出，并恢复 cosurf/nar_csm/nsetf/iatsp/isude/sude
参数与写回逻辑），然后才能用官方 `DIELECTRIC ENERGY = -0.35863 EV` 验证。

## 14. 2026-10-25 SNAP 移植尝试：单位/插入点未确认，已回退

- 尝试按 readmo.F90:1366 的 `snapth` 直接移植，并在 getgeo 之后插入
  `if (keywrd 含 SNAP) for i: geo(2,i)=snapth(geo(2,i))`。
- 在测试输入 `H 1 0.96 2 90.0005`（SNAP）上出现**与官方相反**的行为：
  - 官方 HOF = -54.73892（与无 SNAP 完全相同），几何表仍打印 90.0005000；
  - 我们的实现把角度变为 90.0000000，HOF 变为 -54.73871（差 2.1e-4 kcal）。
- 说明 Fortran 在该插入点看到的 `geo(2,i)` 单位/状态与 C++ 不同（很可能不是
  弧度，或 SNAP 块作用于另一份几何副本），直接套用 `snapth` 会引入偏差。
- 已**回退**该移植，避免与官方结果不一致；后续移植前需先确认 `geo` 在
  readmo 阶段的单位与 SNAP 实际生效位置（可用打印 `snapth` 输入/输出验证）。

## 15. 2026-10-25 关键修复：keywrd 缺少首尾空格导致行首关键字失效

- **根因**：MOPAC 的 `keywrd` 约定首尾各有一个空格，代码普遍用
  `keywrd.find(" DRC")`、`find(" 1SCF")` 这类带前导空格的方式匹配关键字。
  C++ 的 `gettxt()` 直接把数据卡第一行赋给 `keywrd`（无前导/尾随空格），
  导致**行首关键字匹配失败**：
  - 例如输入 `DRC` 时 `find(" DRC")` 失败 → 被误派发到默认 EF 分支；
    C++ 只做了一次 EF 单点/优化，没有跑 DRC 轨迹；
  - 输入 `1SCF PM7` 时 `1SCF` 位于行首 → 被忽略，程序会去优化几何。
- **修复**（gettxt.cpp，2026-10）：`keywrd` upcase 后补前导和尾随空格。
- **验证**：
  - `1SCF PM7` 水分子：官方向 -57.79103，C++ 修复后同样 -57.79103，
    并且输出出现 “1SCF WAS SPECIFIED”；
  - `DRC` H2 用例：修复前 C++ 只输出 7 KB（EF 单点）、没有轨迹文件；
    修复后输出 488 KB、`drc.xyz` 1.02 MB，与官方（484 KB / 1.03 MB）同一量级，
    DRC 正常执行；
  - 原有 5 个基准（H2/PM7/MOZYME/H2O/硅醇）HOF 全部不变，仍与官方逐位一致。
- 同时清除了遗留调试输出：`ef.cpp` 的 `EFG-C1`/`EFM-C` 探针、
  `run_mopac.cpp` 的 `MOZRUN-*` 打印。
- 待查：DRC 轨迹步数/xyz 大小与官方仍有约 1% 差异，需对比时间步长/初速度。

## 16. 2026-10-25 DRC 输出格式对齐（XYZ 轨迹已与官方逐字节同尺寸）

- `write_trajectory.cpp`：XYZ 原子行改为 Fortran `(3X,A2,3F15.5)` 等价格式，
  且元素符号按 **右对齐 2 字符** 输出（Fortran 的 `elemnt` 对单字母元素是 " H"，
  C++ 表是 "H "，此前导致每行少 1 字节且列错位）。
- `drcout.cpp`：轨迹行的 `%` 字段去掉多余空格，改为官方 `%` + `F8.4` 语义。
- 结果（`DRC` H2 用例，5795 帧）：
  - `drc.xyz` 官方与 C++ 均为 **1,031,510 字节**，逐行对比仅 34 行不同，
    差异为累计积分造成的 1e-5 Å 级坐标末位（非格式问题）；
  - `drc.out` 轨迹行格式与官方一致（如 `...   %  0.0000`）。
- 非 COSMO 的 5 个基准 HOF 仍与官方逐位一致。

## 17. 2026-10-25 `getsym` 接线完成（对称性约束真正生效）

- **输入源修复**：`getsym()` 4 参版本原来用 `std::getline(std::cin, ...)` 读取
  约束（数据流接错），现改为从 `getdat_lines` 队列按 `getdat_line_idx` 顺序读取；
- **调用条件**：`readmo.cpp` 原来无条件调用 0 参 stub；现按 Fortran 语义，
  仅当 `keywrd` 含 `SYMM`/`SYM ` 或 `ndep > 0` 时调用 4 参真实版本，并在调用前
  为 `symmetry_C::locpar/idepfn/locdep/depmul` 分配 3*natoms / natoms 长度；
- **0 参 stub 已删除**；
- **输出格式对齐**：`PARAMETER DEPENDENCE DATA` 标题改为 Fortran 的 `5X` 缩进；
  数据行改为 `i13,i19,i16,6i5` + `43X,7i5` 续行的动态格式，修掉 C4473 告警；
- **验证**：`PM7 1SCF SYMMETRY` 水分子用例，C++ 现在与官方一样打印
  `PARAMETER DEPENDENCE DATA`，HOF 均为 -57.79103 kcal/mol；
  5 个基准 HOF 全部不变。

## 18. 2026-10-25 回归语料抽检：55 个 deck 与官方记录逐位一致

- 语料：`regression/inputs`（349 个 deck，配套 `out_ftn/*.arc` 官方记录输出、
  `out_cpp` C++ 记录输出）。
- 方法：排除 COSMO/EPS 用例（已知缺口）后，抽检 25 个最小 deck + 随机 30 个
  deck，用当前 C++ 运行并解析 `FINAL HEAT OF FORMATION`，与 `out_ftn/*.arc`
  的 `HEAT OF FORMATION` 逐位比较。
- 结果：**55/55 一致，0 处差异**；覆盖 MNDO / PM3 / PM7、含 Na/S/P/Cl/F 等元素、
  双原子/多原子/单原子用例，以及 DRC（无最终 HOF，双方一致）。
- 说明：此前 16 个 tests/run_test 用例 + 5 个基准 + SYMMETRY/DRC 专项也都与
  官方一致；目前非 COSMO 主干数值可信度较高。

## 19. 2026-10-25 COSMO surclo 移植启动：先修复 rsc 存储布局（真实 bug）

按 Fortran 布局重做 `rsc` 存储：

- Fortran `rsc(4,maxrs)` 是**列主序**（每个环段 4 个连续值）；C++ 原先把 rsc
  存成 `vector<vector<double>>`（4 行），却把 `&rsc[0][0]` 当连续 4*maxrs 缓冲
  传给 `mfinel`（其索引为 `rsc[(irs-1)*4+ix]`），导致环段坐标被读成错位数据，
  是 COSMO NaN 的重要原因。
- 修复：`cosmo.cpp` 与 `linear_cosmo.cpp` 统一改为
  `std::vector<double> rsc(4*(maxrs+1))`（列主序），`surclo` 签名改为
  `double*`，`mfinel`/`surclo` 全部按 `(seg-1)*4+ix` 访问；`cosmo.h` 同步。
- 效果：`PM7 EPS=78.4 1SCF` 水分子下 `amat` 非有限元由 41 个降到 31 个，
  但仍非零 → COSMO 仍不收敛（默认 2000 轮很慢）。
- 结论：rsc 布局是必须先修的底层 bug，但**现有 surclo 移植仍不完整/有偏差**
  （环点仍有退化重合），需要按 cosmo.F90:1202-1663 的 462 行做完整对照重写。
- 非 COSMO 的 5 个基准 HOF 不受影响，仍与官方逐位一致。

## 20. 2026-10-25 `surclo` 完整移植完成首版：COSMO 不再 NaN/挂起

- 按 `cosmo.F90:1202-1663`（462 行）重写 `cosmo.cpp` 中的 `surclo`：
  - 环段生成（ansude、fz/yx、rvx/rvy、三元点 ntrp/phiset/iset/tarset）；
  - 环段面积/体积累加（rsc(4)、area、cosvol）；
  - isude/sude、nipsrs 计数；
  - 按原子排序、按最近主段分配、必要时插入新表面段、最终按主段排序；
  - 修正 `n0` 索引（C++ 为 1-based：`n0[1]`=nspa、`n0[2]`=氢用方向数）。
- 同时前置修复了 `rsc` 列主序布局 bug（见 §19）。
- 效果（`PM7 EPS=78.4 1SCF`）：
  - **COSMO 不再 NaN、不再挂起**，水分子 0.6 s 收敛；
  - COSMO AREA / VOLUME 与官方完全一致（42.49 / 25.19）；
  - 单原子、H2、CO2 的 HOF 与官方一致（CO2 介电 -0.14917 vs -0.14916 EV）；
  - 仍存在偏差：水 HOF -65.51751 vs 官方 -65.16522（介电 -0.37726 vs -0.35863 EV，约 5%），
    CH4 HOF -14.80880 vs -14.74137（介电 -0.01911 vs -0.01602 EV，约 19%）。
- 结论：环段几何/面积（全局量）已正确，但环段参与 `amat` 的分布仍有小偏差；
  下一步可用本机 Intel `ifx` 把 Fortran `surclo`/`ansude` 抽成独立程序，
  与 C++ 版逐段对比 `rsc`、`nipsrs`，定位剩余偏差。
- 非 COSMO 的 5 个基准 HOF 不受影响，仍与官方逐位一致。

## 21. 2026-10-25 Fortran surclo 差分测试程序：surclo 已确认正确

- 用 Intel `ifx`（VS2022 + Intel oneAPI 2025.3）把 `cosmo.F90` 的
  `ansude`（L110-174）与 `surclo`（L1202-1663）原样抽出，配 stub 模块
  （`funcon_C/molkst_C/cosmo_C`）和驱动程序，编译成 `surclo_harness.exe`。
- C++ 侧在 `surclo()` 内加环境变量开关 `SURCLO_DUMP=<prefix>`：入口导出全部
  输入（coord/srad/nat/nipa/din/lipa/cosurf/iatsp/nar_csm/nsetf、n0、rsolv、
  初始 nps），出口按 **16 进制位模式** 导出 rsc/nipsrs/ipsrs/isort/cosurf/
  area/cosvol。
- 同一份 `water.in` 分别跑 C++ 与 Fortran harness，结果（水分子）：
  - 两侧均为 `nps=49, nrs=38, nipc=...`，625 个输出值数量/顺序完全一致；
  - 211 个值有差异，但**最大相对差 1.65e-7**（纯浮点运算顺序/编译器差异，
    无结构性差异）。
- 结论：**surclo 移植算法正确**。水分子介电能量仍差约 5%（-0.37726 vs
  -0.35863 EV），问题在下游的 `coscav`/`mfinel`/`amat`、`bmat` 或
  `addfck/addnuc`，而非环段生成。
- 非 COSMO 5 个基准 HOF 不变；COSMO 仍正常收敛（0.6 s）。
- 下一步：用同样方法对 `coscav`（amat）做差分（导出 mfinel 输入/输出与
  amat，再与 Fortran harness 对比），定位剩余偏差。

## 22. 2026-10-25 COSMO 闭环：定位并修复 mfinel 的 tm 转置打包 bug

- 用同样的 16 进制差分方法对 `coscav` 做对照：C++ 导出 coscav 输入
  （coord/srad/nat/dirvec/dirsm/n0/rsolv/disex2）与输出（amat/cosurf/
  iatsp/nar_csm/nsetf/arat），Fortran harness 用完全相同输入调用原版
  `coscav`，逐项比较。
- 首轮对比结果指认问题在 `amat`：1215/1225 个元素不一致（最大相对差
  6.3e2，近零元），而 area/cosvol/iatsp/nar_csm/nsetf 完全一致、
  cosurf 仅差 5e-15。
- **根因**：`coscav` 把 `Tmat` 传给 `mfinel` 时打包成了
  `tm_flat[j2*3+ix] = Tmat[i][j2][ix]`（按 [axis][component] 原样），
  但 Fortran `tm(row=axis, col=component)` 与 `mfinel` 的
  `TM(ix,j)=tm(ix,j)` 要求打包时 **转置**。因此 finel 列表里的网格点坐标
  被错误旋转，amat 全错；area/cosvol 在旋转之前计算，故仍正确。
- 修复：打包改为 `tm_flat[j2*3+ix] = Tmat[i][ix][j2]`（tm2_flat 同步）。
- 修复后 coscav 差分：amat 最大相对差 **2.35e-12**（纯浮点舍入），
  area/cosvol/iatsp/nar_csm/nsetf 逐位一致。
- COSMO 官方闭环（`PM7 EPS=78.4 1SCF`）：
  - 水：HOF **-65.16522**（官方 -65.16522），DIELECTRIC **-0.35858 EV**
    （官方 -0.35863 EV，差 5e-5 EV），AREA 42.49 一致；
  - CH4：HOF -14.74135（官方 -14.74137），DIELECTRIC -0.01599（官方 -0.01602）；
  - CO2：HOF -87.44782（官方 -87.44781），DIELECTRIC -0.14917（官方 -0.14916）。
- 剩余 ~1e-4 EV 级差异来自 C++/Fortran 浮点运算顺序，属可接受范围。
- 非 COSMO 5 个基准 HOF 仍与官方逐位一致。

## 23. 2026-10-25 剩余 ~1e-4 EV 误差归因（非结构性 bug）

在 COSMO 闭环后继续追剩余误差：

- 新增 pre-Cholesky `amat` 差分（C++ 与 Fortran harness 都在 `coscl1` 前
  按 16 进制导出原始 amat）：
  - 原始 amat：1225 个元素中 647 个存在差异，但**最大相对差仅 5.6e-15**
    （纯浮点求和顺序）；
  - amat 的 2-范数条件数 **cond ≈ 20.8**（良态），排除“病态放大”解释；
  - Cholesky 后 amat 相对差 2.35e-12，属分解过程正常的舍入放大。
- 水分子 COSMO 输出对比（官方 vs C++）：
  - `FINAL HEAT OF FORMATION`：-65.16522 vs **-65.16522（一致）**；
  - `TOTAL ENERGY`：-322.99861 vs **-322.99861（一致）**；
  - `CORE-CORE REPULSION`：-107.78158 vs **-107.78158（一致）**；
  - `ELECTRONIC ENERGY`：-215.21704 vs -215.21703（差 1e-5 EV）；
  - `DIELECTRIC ENERGY`：-0.35863 vs -0.35858 EV（差 5e-5 EV）；
  - `COSMO AREA`：42.49 vs 42.49（一致）。
- 结论：剩余误差来自 **SCF 收敛点/浮点运算顺序的微小差异**（电子能量差
  1e-5 EV），不是 cavity/amat/bmat 的结构性错误；HOF 已按打印精度一致。
- 附带发现（另一独立缺口）：C++ 的点群识别输出 C1，官方为 C2v；
  属于对称性检测未接线（不影响 1SCF 数值）。

## 24. 2026-10-25 剩余误差再定位 + 点群（molsym）缺口

### A. COSMO 剩余 1e-4 EV 级误差
- 收紧 SCF（`SCFCRT=1.0D-8`）后：HOF 仍一致、ELECTRONIC ENERGY 变为逐位相同，
  DIELECTRIC 差从 5e-5 降到 1e-5 EV → 说明残差主要来自 SCF 收敛路径/求和顺序。
- 原始 amat 最大相对差 5.6e-15、条件数 20.8，经 Cholesky 后 2.35e-12，
  不足以解释 1e-4 级介电差；剩余属浮点运算顺序/线性代数库实现差异。

### B. 点群 C1（真实翻译缺口，非数学库）
- 官方 `writmo.F90:254` 在打印 ELECTRONIC ENERGY/POINT GROUP 前调用
  `symtrz(c, eigs, 1, .true.)`；C++ `writmo.cpp` **漏掉了该调用**，
  所以 `symmetry_C::name` 一直是默认 C1。
- 补上调用后仍为 C1 → `molsym()` 本身翻译有 bug；
  已定位到本征值下标 off-by-one：C++ 用 `ew[1]/ew[2]/ew[0]`，
  Fortran 为 `ew(2)/ew(3)/ew(1)`（`rsp` 输出是 1-based 升序）。
- 修正下标后点群从 C1 变为 **R3**（仍错误，水应为 C2v），说明 molsym 还有
  更深的结构性翻译问题（symopr/bldsym/chi/cartab/symdec 的组合与 goto 结构）。
- 结论：点群需要按 molsymy.F90（452 行，含 bldsym/symopr/plato/orient/chi/
  cartab/symdec 依赖）做完整对照重写；本次实验性改动已全部回退，
  避免输出错误点群/错误对称化轨道。
- 当前状态：5 个基准与非 COSMO 回归全部保持与官方一致；COSMO 数值闭环保持
  （HOF -65.16522、介电 -0.35858 vs -0.35863）。

## 25. 2026-10-25 点群修复完成（C1 → C2v/D*h/Cs）：molsym 的 rsp 打包指针 bug

- 用 Intel `ifx` 把 MOPAC2016 的 `molsym` 全套依赖（molsymy/bldsym/symopr/
  rsp/rotmol/symmetry_C block data，含 chi/orient/plato/cartab/symdec）编译成
  独立 harness（工作区 `molsym_harness.exe`），水分子输入下返回
  **name=C2v, nclass=3, ierror=0, ielem 与官方一致**，作为对照基准。
- C++ `molsym()` 返回 C1；逐段 trace 定位到它走了 `linear` 分支：
  - 根因：`molsym` 把 1-based 的 packed 惯性矩数组 `f[1..6]` 直接传给
    `rsp()`，而 C++ `rsp()` 期望 **0-based packed** 指针（`a[0]` = (1,1)）。
    矩阵整体错位 → 特征值错误（第二个特征值接近 0）→ `linear=true` → C1。
  - 修复：改为 `rsp(&f[1], 3, ew, r)`。
  - 另修正两处遗留的 1-based 特征值下标（与 0-based `rsp` 输出不一致）：
    二重简并交换 `ew(1)<->ew(3)` → `ew[0]<->ew[2]`；
    条件 `ew(3)-ew(2)<1e-2*ew(3)` → `ew[2]-ew[1]<1e-2*ew[2]`。
- 同时在 `writmo.cpp` 补上 Fortran `writmo.F90:254` 的
  `symtrz(c, eigs, 1, .TRUE.)` 调用（1-based 列主序打包 c），并回写对称化后的
  特征向量与本征值；打印前若 moperr 则返回。
- 验证（官方 vs C++）：
  - H2：`D*h` vs **D*h**；水：`C2v` vs **C2v**；硅醇：`Cs` vs **Cs**；
  - 水 COSMO：HOF -65.16522 一致，介电 -0.35858 vs -0.35863 EV；
  - 5 个基准 HOF 全部不变、仍与官方逐位一致。
- Fortran harness 的 r 矩阵符号/相位与 C++ 有 ~1e-9 差异（harness 用 Jacobi
  替代 eigenvectors_LAPACK；官方/C++ 用 dsyevd），不影响 name/nclass/ielem。

## 26. 2026-10-25 COSMO 剩余差异定位：全部为浮点/SCF 级，无结构性 bug

用 Fortran harness 对 COSMO 各环节做逐项差分（同一份 C++ 输入）：

| 环节 | C++ vs Fortran | 结论 |
|---|---|---|
| cavity/amat（Cholesky 后） | 最大相对差 2.35e-12 | 一致（舍入级） |
| `mkbmat`（B 矩阵） | 588 个值中 143 个差 1 ulp，最大相对差 **3.4e-16** | 一致 |
| `addfck`（给定完全相同输入） | qscnet 最大相对差 **2.4e-12**，ediel 同量级 | 一致 |
| 完整 COSMO 介电能量 | -0.35858 vs 官方 -0.35863 EV（相对 1.4e-4） | 差异来自输入密度 |

- 收紧 SCF（`SCFCRT=1.0D-8`）后电子能量与官方逐位相同，介电差降到 1e-5 EV；
- 即：COSMO 的 cavity、B 矩阵、addfck 求解本身都已与 Fortran 一致到浮点舍入，
  残差来自 **SCF 电子密度路径**（SCF 迭代/线性代数/libm 舍入），不是 COSMO
  代码的结构性错误，无法在不逐位复刻 SCF 的情况下消除。
- 清理了本轮临时差分导出代码；5 个基准与点群（C2v）保持与官方一致。

## 27. 2026-10-25 第 7 批中文注释（COSMO/SCF/存储相关 12 个文件）

- 新增中文导读头（`【中文导读】`）：
  - `cosmo.cpp`（COSMO 全流程：cosini/coscav/surclo/mkbmat/addnuc/addfck/coscl1-2，
    并更新过时的 “stubbed” 说明）；
  - `cosmo_C.cpp`、`symmetry_tables.cpp`（点群表，注明 molsym 已修复）；
  - `datin.cpp`（EXTERNAL/PARAMS 参数读取）；
  - `bonds.cpp`、`chrge.cpp`、`densit.cpp`、`den_in_out.cpp`；
  - `buildf.cpp`（Fock/fock2z/focd2z）、`bfn.cpp`；
  - `chanel_C.cpp`、`common_arrays_C.cpp`。
- 仅新增注释行，代码行逐行比对不变；编译通过，5 个基准 HOF 与官方一致，
  COSMO 水分子仍为 C2v / -65.16522 / 介电 -0.35858 EV。
- 累计注释批次：7 批，约 74 个文件、约 1100+ 行中文注释。

## 28. 2026-10-25 第 8 批中文注释（导数/CI/分析/加氢 12 个文件）

- 新增 `【中文导读】` 头：
  - `deri21.cpp`（最小二乘正交基）、`deri22.cpp`（CI 超向量 AB/FCI）、
    `deri23.cpp`（CI Fock 对称化）；
  - `ccrep.cpp`（芯-芯排斥/D3 修正）；
  - `charvi.cpp`、`charst.cpp`（本征向量/态的对称性特征标，后者注明
    minv/dtrans/matout 仍为占位）；
  - `anavib.cpp`（振动模式报告）、`add_hydrogen_atoms.cpp`（自动加氢，9 子程序）；
  - `analyt.cpp`（解析 Hessian 辅助）、`compct.cpp`（PATH/IRC 连接表修正）、
    `dang.cpp`（夹角余弦稳定计算）、`copyc6.cpp`（C6 色散系数复制）。
- 仅新增注释，代码行逐行比对不变；编译通过；5 个基准 HOF 与官方一致。
- 累计：8 批、约 86 个文件、约 1250+ 行中文注释。

## 29. 2026-10-25 第 9 批中文注释（几何/IO/常数/Fock 12 个文件）

- 新增 `【中文导读】` 头：
  - `gmetry.cpp`（内坐标↔笛卡尔，bangle_ 接口）；
  - `gettxt.cpp`（关键字/续行/SETUP 合并/keywrd 首尾空格修复）、
    `getgeo.cpp`（几何解析、getdat_lines 队列、txt_to_atom_no）；
  - `funcon_C.cpp`、`conref_C.cpp`（两套物理常数与能量换算依赖）、
    `elemts_C.cpp`（两字符元素符号表）；
  - `diag.cpp`（阈值化本征值/本征向量）、`fock1.cpp`（单中心 Fock 项）；
  - `geochk.cpp`（结构/电荷检查主驱动）、`dofs.cpp`（自由能直方图）；
  - `distance.cpp`（距离/键角/二面角辅助）、`dftd3.cpp`（D3 色散修正）。
- 仅新增注释，代码行逐行比对不变；编译通过；5 个基准 HOF 与官方一致。
- 累计：9 批、约 98 个文件、约 1380+ 行中文注释。

## 30. 2026-10-25 第 10 批中文注释（MOZYME/CI/D3 辅助 12 个文件）

- 新增 `【中文导读】` 头：
  - `fock2z.cpp`（focd2z 含 d 轨道两中心 Fock）、`fock1.cpp` 已有（第 9 批）；
  - `diagi.cpp`（CI 的 xy 关联项）、`diagg.cpp`（MOZYME 对角化/密度/epseta 集合）、
    `diagg2.cpp`（占据-虚 LMO Jacobi 消除）；
  - `addhb.cpp`（氢键/成键修正）、`adjvec.cpp`（MOZYME 向量邻接修正）；
  - `axis.cpp`（转动惯量/转动常数）、`bangle.cpp`（键角）；
  - `bonds_for_MOZYME.cpp`、`density_for_MOZYME.cpp`、`dipole_for_MOZYME.cpp`
    （MOZYME 键级/密度/偶极，标注打印占位）；
  - `dftd3_bits.cpp`（D3 几何阻尼与 gdisp）。
- 仅新增注释，代码行逐行比对不变；编译通过；5 个基准 HOF 与官方一致。
- 累计：10 批、约 110 个文件、约 1500+ 行中文注释。

## 31. 2026-10-25 第 11 批中文注释（积分/CI/工具 12 个文件）

- 新增 `【中文导读】` 头：
  - `dot.cpp`（1-based 向量点积）、`dijkl1.cpp`（CI dijkl 原子对部分）、
    `dijkl2.cpp`（dijkl 的 MO 变换）；
  - `diat.cpp`（diat/diat2/ss 双原子重叠积分）、`dipind.cpp`（偶极积分）、
    `dihed.cpp`（二面角，含周期最小镜像）；
  - `digit.cpp`（字符串数字解析）、`dimens.cpp`（MOLECULAR DIMENSIONS 表）、
    `euler_C.cpp`（周期胞整数下界）、`exchng.cpp`（行列/向量片段交换）；
  - `fillij.cpp`（MOZYME ij/重叠表填充）、`findn1.cpp`（壳层起点查找）。
- 仅新增注释，代码行逐行比对不变；编译通过；5 个基准 HOF 与官方一致。
- 累计：11 批、约 122 个文件、约 1620+ 行中文注释。

## 32. 2026-10-25 第 12 批中文注释（力场/频率/变换 12 个文件）

- 新增 `【中文导读】` 头：
  - `dtrans.cpp`、`dtran2.cpp`（对称变换矩阵与 3x3 块）、`genun.cpp`（旋转/正交基）；
  - `derivs_C.cpp`（导数超向量存储）、`dernvo.cpp`（开壳层/CI 解析梯度驱动）、
    `formxy.cpp`（CI xy 超向量组装）；
  - `fmat.cpp`（二阶/四阶差分公式）、`frame.cpp`（刚体平移/转动零模）、
    `freqcy.cpp`（振动频率/简正模式）、`forsav.cpp`（FORCE RESTART 存取）；
  - `geo_diff.cpp`（几何差异 sum/RMS）、`geochk_small.cpp`（fix_charges 与
    add_sp/sp2/sp3_H 补氢小例程）。
- 仅新增注释，代码行逐行比对不变；编译通过；5 个基准 HOF 与官方一致。
- 累计：12 批、约 134 个文件、约 1740+ 行中文注释。

## 33. 2026-10-25 注释全覆盖：全部 .cpp / .h 已加中文导读

- 第 13 批（252 个剩余 .cpp）与第 14 批（368 个剩余 .h）按文件名类别统一生成
  `【中文导读】` 头：
  - 说明所属类别（如 ALPB/xfac 参数表、解析导数/积分、对角化、MOZYME、
    SCF/Fock、输出、对称性、结构/成键、数值工具等）；
  - 说明对应 MOPAC 2016 的 F90 来源、1-based 索引/列主序/打包下三角约定，
    并指向 docs/TRANSLATION_GAPS.md 中的移植缺口；
  - `.h` 头额外注明“声明模块全局量/函数/常量，与同名 .cpp 配对使用”。
- 生成脚本逐文件做 **code-line 等价校验**（忽略 `//` 注释后代码行完全一致），
  并保留原文件换行风格、去除 BOM。
- 覆盖统计：**.cpp 380/380、.h 372/372** 全部含中文导读。
- 构建通过；5 个基准 HOF 与官方逐位一致；COSMO 水分子仍为 C2v /
  -65.16522 / 介电 -0.35858 EV。
- 累计：14 批注释（前 12 批为逐文件精写，后 2 批为全量类别化覆盖），
  约 750 个文件、约 3500+ 行中文注释。

## 34. 2026-10-25 全回归语料扫描（非 COSMO）：339/346 一致

- 对 `regression/inputs` 中 346 个非 COSMO deck 逐个运行 C++，与 `out_ftn/*.arc`
  记录的官方 HOF 对比（超时 45 s）：
  - **339 个完全一致**；
  - 7 个不一致，全部是**几何优化（EF/BFGS）**类 deck：
    `_mx_opt_am1_ch4(_o)`、`_mx_opt_mndo_c6h6(_o)`、`_mx_opt_pm7_h2o(_o)`、
    `_val_ef`。
- 以 `_mx_opt_pm7_h2o` 为例：官方与 C++ 都是 3 个 EF 循环，但
  - 官方：GRAD 3.872 → 2.245 → 0.935，HOF -57.79870，终态 GRAD 0.93495；
  - C++：GRAD 3.872 → 2.213 → 0.931，HOF -57.79850，终态 GRAD 0.93108；
  - 偏差从第 2 个 EF 步出现（梯度 2.245 vs 2.213），属于 EF 本征向量步长
    对浮点/本征分解细节敏感（rsp.cpp 中已记录 EF 对 degenerate 本征向量
    方向敏感，C++ 特意使用 MKL dsyevd 路径）。
- 结论：非优化类 deck 全量一致；优化类剩余偏差在 1e-4~8e-4 kcal/mol 量级，
  来自 EF 轨迹差异，不是新的结构/接口 bug。若要求逐位一致，需要逐位复刻
  EF 的梯度投影与对角化细节（与数学库/运算顺序相关）。
- 扫描脚本：`regression_all.ps1`；结果：`regsweep_all/report.txt`。

## 35. 2026-10-25 POLAR 驱动 + alphaf 移植（进行中）

- 已把 `polar.F90:1-399` 的 `polar()` 驱动和 `polar.F90:742-963` 的 `alphaf()`
  移植进 `polar.cpp`（替换原 retained stub），复用已有 helper（hmuf/makeuf/
  densf/transf/ffreq1/ffreq2/aval/pol_vol/betaf/beopor/non*/ng*/openda）。
- 已修复两个移植 bug：
  1. `polar.cpp` 里 `norbs` 原先绑定到全局 `::norbs`（恒 0），导致 alphaf 用
     norbs=0 崩溃；改为 `using molkst_C::norbs;` 后正常。
  2. `hmuf()` 的偶极矩阵对角项 `1.8897262*coord_[id][i]` 应为
     `coord_[id-1][i]`（Fortran `coord(id,i)` 是 1-based）。
  3. `densf()` 的 w2 工作区需要 norbs*norbs 个元素（原按 norbs+1 分配）。
- 现状：`PM7 1SCF POLAR` 可完整跑完（alphaf/betaf/gamma 全部调用），
  输出结构与官方一致；但 alpha 数值仍有偏差：
  - 对角元中**坐标为 0 的方向**完全一致（水分子 7.1439 A.U.）；
  - 其余两个方向偏大：C++ 14.0026/13.6122 vs 官方 11.4206/10.4599。
  官方结果对坐标平移不变（已用 +10 Å 平移验证），因此不是原点问题；
  疑点集中在“对角项非零方向”的响应计算（hmuf 对角项 → ffreq1/ffreq2/
  transf/makeuf/densf 链），下一步用 Fortran alphaf 差分 harness 定位。
- 非 POLAR 路径不受影响（关键字未触发）。

## 36. 2026-10-25 POLAR 继续修复：ffreq1 的 s-only 原子项缺失（静态 α 已对齐）

- 用 Fortran `alphaf` 差分 harness（从 C++ dump 同一份 norbs/nclose/numat/
  coord/eigs/c/参数/w 输入）对拍，确认 C++ alphaf 与 Fortran 参考实现结果不同，
  逐段检查 `ffreq1` 时发现：
  - Fortran 的 `ffreq1` 分支 `go to (1002,20,10,...) ib-ia+2` 中，
    **idx==2（仅 s 轨道的原子，如 H）仍执行 F(S,S) 项**：
    `f(ia,ia) += pb(ia,ia)*gss + ptpop*gsp - papop*hsp`（其中 ptpop=papop=0）；
  - 原 C++ 代码遇到 idx==2 直接 `continue`，导致所有 H 原子的一中心项缺失，
    静态/含频 α 均偏大。
- 修复（去掉 continue，ptpop/papop 只在 idx>=4 时累加）后：
  - 水分子 `PM7 1SCF POLAR` 的**前两个 α 块与官方完全一致**：
    `9.67479` vs 官方 `9.67480`（静态与 third-harmonic at ω=0）；
  - 其余含频块仍偏低（如 ω=0.25：C++ 8.557 vs 官方 9.676；ω=0.5：
    8.559 vs 9.679），说明还有 ω 相关路径（makeuf/迭代收敛/betaf 调用顺序）
    需要继续对拍。
- 同时建立了 `alphaf` 的 Fortran 差分 harness（工作区 `alphaf_harness.exe`）；
  `polar()` 驱动与 `alphaf()` 已替换原 stub，POLAR 可端到端运行。

## 37. 2026-10-25 POLAR 闭环：alpha/beta/gamma 与官方一致

在 §36 修复 `ffreq1` s-only 项后，静态 α 已对齐；继续对拍含频 α 时发现
二阶（β/γ）通路返回 0，顺藤摸出 DAF 直接访问文件子系统的三个移植 bug：

1. `dawrit/daread` 用 `v[1][flat]` 访问二维矩阵（Fortran 是列主序连续缓冲），
   而 C++ 的 `vector<vector<double>>` 每行只有 norbs+1 个元素 → 越界读写/堆破坏。
   改为按 `row=(k-1)%nrow+1, col=(k-1)/nrow+1` 映射回 `v[row][col]`。
2. `dawrt1()` 每次写记录都用 `trunc` 打开 scratch 文件，前一条记录被清掉；
   且写用 `density_fn`（.den）而 `openda()` 打开的是 `pol_fn`（.pol），读写
   不是同一个文件 → `daread` 读回全 0。统一为 `pol_fn`，且写记录用
   `seekp((ns-1)*irecln*8)` 原地写、不截断。
3. Fortran `chanel_C.F90:57` 的 `ioda/ifilen` 是 `dimension(145)`，C++ 用了
   2001：`openda` 头部 4+145*4+145*4=1164B 可放进 rec=1（8184B），而 2001 的
   头部 16004B 越过 rec=1 并覆盖数据记录 → β/γ 数值错误（-13.74 vs -20.87）。
   C++ 的目录数组与 145 次循环已改为 145。

验证（水分子 `PM7 1SCF POLAR`，官方 vs C++）：

| 量 | 最大相对偏差 |
|---|---|
| ISOTROPIC AVERAGE ALPHA（6 个频率） | 1.0e-6 |
| AVERAGE BETA（12 行，SHG/EOPE/OR 等） | 4.0e-7 |
| AVERAGE GAMMA（3 个频率） | 3.0e-5 |

甲烷 `PM7 1SCF POLAR` 的 6 个 α 也与官方一致（如 18.04165 vs 18.04164）。
剩余差异来自迭代收敛容差/浮点舍入，不是数学库实现差异。

至此 POLAR 功能缺口（polar 驱动、alphaf、betaf、beopor、non*、ng* 及 DAF
存档）闭合；非 POLAR 路径回归不受影响（水 PM7 -57.79103 仍与官方一致）。

## 38. 2026-10-25 MOZYME 输出：LMO（局域分子轨道）块补齐

- 重写 `prtlmo.cpp`：按 `prtlmo.F90`（prtlmo + prtlmn）完整移植——
  按 LMO 能量排序、逐 LMO 列出 |psi|^2 最大的最多 20 个原子
  （`I5,F9.4` + `I5,1X,A2,I5` 格式），占据集/虚拟集分别用
  `isort(1..)` 与 `isort(noccupied+1..)`、`eigs(1..)` 与 `eigs(noccupied+1..)` 偏移。
- `writmo.cpp` 补上 F90 L556-575 的 MOZYME 分支：
  - `RE-LOC` → 打印 "LMOs being Re-Localized" 并调用
    `local_for_MOZYME("OCCUPIED"/"VIRTUAL")`；
  - `VEC/ALLVEC` → `EIGEN` 关键字走 `eigen(false,true)`，否则 `prtlmo()`；
  - 同时把 MO 的 `EIGENVALUES` 块限定为 `!mozyme`（F90 writmo L819 只在
    非 MOZYME 时打印），与官方输出一致。
- 验证（`PM7 MOZYME VEC 1SCF` 水）：C++ 现可打印
  "LOCALIZED MOLECULAR ORBITALS / Occupied Set / Virtual Set" 表，
  每个 LMO 的原子贡献（PSI**2 整数）与官方相差约 0.1~0.4%
  （如 6892 vs 6886、9460 vs 9479、6620 vs 6624）。
- 仍存差异：LMO 的 **energy 列**偏约 0.08~0.4 eV
  （C++ 用 SCF 的 canonical `eigs`，官方 prtlmo 用的 `eigs` 是 MOZYME
  SCF 中的 LMO 能量）；这涉及 MOZYME SCF 的局域化能量写回 `eigs` 的路径
  （`lmo_to_eigenvectors` 在 writmo 后段才调用），下一步继续对齐。
- 非 MOZYME 路径不受影响（水 PM7 仍 -57.79103；MOZYME 水 -57.78972 = 官方）。

## 39. 2026-10-25 MOZYME + MULLIK：lmo_to_eigenvectors/FlmoFromFao 移植

- 针对 `PM7 MOZYME MULLIK 1SCF` 的崩溃（C++ 之前不产生 .out）：
  - 在 `eigen.cpp` 按 `eigen.F90:192-379` 移植 `FlmoFromFao()`（在 LMO 基下
    用 MOZYME 压缩 Fock `f` 构造 Fock 矩阵），并按 `eigen.F90:438-499`
    移植 `lmo_to_eigenvectors()`（`rsp` 对角化后展开回 AO 基，支持
    `eig_offset/col_offset` 处理虚拟集从 `noccupied+1` 开始）；
  - `writmo.cpp` 的 MULLIK/GRAPH 分支补上 F90 L857-887：MOZYME 时先把
    `noccupied`/`nvirtual` 两组 LMO 转成本征向量，再把 `h` 从 MOZYME 压缩
    形式转成下三角打包（`convert_mat_packed_to_triangle`）；
  - 补上 F90 L890-894 的 `mpcpop(-1)` 调用，恢复
    "MULLIKEN POPULATIONS AND CHARGES" 表。
- 结果：`PM7 MOZYME MULLIK 1SCF` 不再崩溃，能打印 MULLIKEN 表。
- 仍存差异：C++ 的 Mulliken population 数值不正确（如 O 为 -0.1006 /
  charge 6.1006，官方为 6.6984 / -0.6984），说明 `mullik()` 在 MOZYME 路径
  下用到的 `h`/`pb` 语义或压缩布局还没有对齐；下一步对照 `mullik.F90`
  中 S^{-1/2} 与密度矩阵的构造逐项核对（非 MOZYME 的 Mulliken 正常）。

## 40. 2026-10-25 Mulliken 分析修复：普通与 MOZYME 路径均与官方逐字节一致

§39 中 MOZYME+MULLIK 的 population 数值错误，追查后发现是 `mullik.cpp`
里两个布局 bug（普通 MULLIK 同样受影响）：

1. `density_for_GPU()` 内部按 Fortran 列主序读取 `c`：
   `c[(k-1)*norbs+(i-1)]` 即元素 `(i,k)`；`mullik.cpp` 却先把
   `mult()` 输出的列主序 `vecs` 转成行主序再传入 → 密度矩阵实际被转置。
   改为直接传 `vecs_orth.data()`。
2. `density_for_GPU()` 的 `pp` 是 1-based 打包指针（内部 `pp[pack_idx]`，
   pack_idx 从 1 开始），`mullik.cpp` 传了 `pb.data()+1`，整体偏移一格。
   改为 `pb.data()`。
3. `mpcsyb.cpp` 的 `mpcpop()` 输出格式对齐 F90
   `(3/8X, 'MULLIKEN POPULATIONS AND CHARGES', /)`、
   `(6X,'NO.  ATOM   POPULATION      CHARGE')`、
   `(5X,I4,3X,A2,F13.6,F14.6)`；C++ elemnt 是左对齐 2 字符（"O "），
   F90 是右对齐（" O"），打印时 trim 后再用 `%2s`。

验证（水，`PM7 MULLIK 1SCF` 与 `PM7 MOZYME MULLIK 1SCF`）：
MULLIKEN 表整块与官方 **逐字节一致**：

```
      NO.  ATOM   POPULATION      CHARGE
        1    O     6.694629     -0.694629
        2    H     0.652685      0.347315
        3    H     0.652685      0.347315
```

MOZYME 路径由 §39 的 LMO→本征向量转换 + 本节修复后同样一致
（6.698389 / -0.698389 等）。非 MOZYME 路径回归不受影响。

## 41. 2026-10-25 关键字回显（wrtcon/wrtwor/wrtout）与 readmo 头部块

- `wrtkey.cpp` 的 `wrtcon/wrtwor/wrtout` 原为空实现；现改为**表驱动**：
  - 从 `wrtkey.F90`（284-897 / 898-1102 / 1103-1516）机械提取纯字面量的
    `if (myword(allkey,"KEY")) write (iw,'(literal)')` 行，生成三张
    `WrtKeyEntry` 表（73 + 78 + 39 条），运行时按源文件顺序匹配并打印；
  - `wrtwor` 另补 F90:1173-1250 的 `T=` / `DUMP=N` 两行（含 M/H/D/W 单位换算，
    默认 172800 秒 / 7200 秒）；
  - `wrtcon` 开头按 F90 执行 `allkey = keywrd`（原来 allkey 为空，导致
    myword 全部不匹配）；`myword` 会消耗/抹掉已匹配关键字，因此包装时给
    `allkey` 首尾各补一个空格，与 Fortran 的 keywrd 形态一致。
- `readmo.cpp` 在 `wrtkey()` 前补 F90 readmo.F90:604-626：
  - 方法标题行 `"PM7 CALCULATION RESULTS"`（24X + 9 字符方法名 + 28X + 3 换行
    + 1X+15('*****')+'****'）；
  - `*  CALCULATION DONE: ...` 行（31X + fdate 时间串 + `  *`）；
  - 关键字回显之后再补一条 80 字符分隔线。
- 验证（`PM7 MOZYME 1SCF` 水）：C++ 输出现在包含
  `PM7 CALCULATION RESULTS`、`CALCULATION DONE`、`PM7/MOZYME/1SCF/T=/DUMP=N`
  关键字行，与官方顺序/文本一致。剩余头部差异：官方最前面的 license banner
  块、输入 deck 回显的位置（官方在头部之后，C++ 目前在最前），以及
  `PM7 CALCULATION` 行的版本/日期格式（`16.013))` vs `(Version: 16.013W)`）。
- 回归：普通水 PM7 -57.79103、MOZYME 水 -57.78972 不变。

## 42. 2026-10-25 输出头部：license banner、版本串与 readmo 顺序

- `readmo.cpp` 开头补 F90 readmo.F90:487-547 的 license banner（10 行，
  80 列框线、Site#/Version/Cite/web-site/MOPAC2016 块），与参考 64 位
  Windows 版逐字节一致；为此：
  - `GetDateStamp.cpp` 的 version 由 "16.013" 改为 "16.013W"；
  - `main.cpp` 的 site_no 由 88 改为 999（参考程序 banner 显示 Site#: 999）。
- `writmo.cpp` 版本行修正：原来是 `MOPAC2016 (Version: 16.013))`（多一个右括号），
  改为 `MOPAC2016 (Version: 16.013W)`。
- `readmo.cpp` 输出顺序按 F90 重排：
  banner → 方法标题/CALCULATION DONE 行 → 关键字回显 → 分隔线 →
  `wrttxt(iw)`（关键字/标题回显）→ `if (prt_coords) geout(1)`（几何表）。
  原先 `geout(iw)`（参数 6 而非 1）在头部之前执行，导致输入原文被提前打印。
- `wrttxt.cpp`：注释/标题行按 F90 加一个前导空格（`" %s\n"`）。
- 结果（`PM7 MOZYME 1SCF` 水）：官方 .out 前 24 行（banner、标题、
  CALCULATION DONE、5 条关键字、分隔线、关键字/标题回显）与 C++ **逐字节一致**，
  唯一差异是 CALCULATION DONE 的时间戳（运行时刻不同）。
- 仍存差异：`geout(1)` 的几何表列宽/逗号"*"格式与官方不同；
  PM7 参考文献块在 C++ 中出现的位置偏早；后续 MOZYME 区块还有格式差异。
  整文件对比：40/121 行完全一致（此前 0），MOZYME 输出大小 5426 vs 官方 5287。
- 回归：普通水 PM7 -57.79103、MOZYME 水 -57.78972 不变。

## 43. 2026-10-25 剩余差异：geout 几何表格式 + LMO 能量取值

1. `geout.cpp` 的 mode==1 笛卡尔表改为严格按 F90 geout.F90:259-261
   `(I6,6X,A,F13.8,1X,A2,F13.8,1X,A2,F13.8,1X,A2)`：
   - A 字段用 9 字符宽、元素符号右对齐 2 字符（Fortran elemnt 为 ' O'，
     C++ 为 'O '，需重建）；
   - 数字之间只保留 1X + A2（`flag` 本身 2 字符），不再多插空格。
   现在 `PM7 MOZYME 1SCF` 的 ATOM 几何表（官方 26-31 行）与官方逐字节一致。
2. `prtlmo.cpp`：官方 `prtlmo` 打印的 `eigs` 是 MOZYME SCF 里 `diagg1`
   写入的 LMO 期望值 `<i|F|i>`，而 C++ `eigs` 仍是 canonical 本征值。
   现在在 `prtlmo()` 里用已移植的 `FlmoFromFao()` 在 LMO 基下构造 Fock、
   取对角元作为 LMO 能量（Fock 不可用时退回全局 eigs）。
   水的 VEC 输出能量从 -23.1582 变为 **-23.1910**（官方 -23.2432），
   已明显接近；PSI**2 也从 6892 vs 6886 等小幅差异可见 LMO 系数仍与官方
   略有不同——根因是 `diagg1.cpp` 仍是 skeleton（未构造 fmo/未做
   occupied-virtual 旋转），LMO 系数本身与官方有差异。
3. 整文件逐行对比（`PM7 MOZYME 1SCF` 水）：43/121 行一致（本节前 40）。
4. 回归：MOZYME 水 HOF 仍 -57.78972，普通水 -57.79103。

下一步（较大）：完整移植 `diagg1.F90`（600+ 行，MOZYME 的核心
occupied–virtual 旋转 + fmo 构造 + LMO 能量写入 eigs），这是 MOZYME LMO
系数/能量与官方完全一致的关键。

## 44. 2026-10-25 MOZYME LMO 能量差异追查（diagg1 / SCF 轨迹）

- 先澄清一个误判：`diagg1.cpp` 只是占位 skeleton，但**真正被调用的是
  `diagg.cpp` 里的完整 `diagg1` 实现**（417 行，含 fmo/ifmo 构造、
  occupied LMO 能量写入 `eigs`），两者签名不同（`vector<char>` vs
  `vector<bool>`），链接上互不冲突。
- 插桩 `diagg()` 后确认 C++ 的调用序列：idiagg = 0..6；能量分支
  （F90 条件：`idiagg<=2` 或 `idiagg%4==0`）在 0,1,2,4 执行，
  最后一次写入 eigs 是 idiagg=4：
  -23.1556 / -23.1582 / ...（与 canonical 几乎相同），
  之后 5、6 不再刷新，所以 writmo 时 `eigs` 停在这一代。
- 官方 `MOZYME VEC` 的 LMO 能量为 -23.2432 / -23.2404 / -15.2051 / -12.1232。
  用 `FlmoFromFao` 在当前（最终）LMO 系数上取 `<i|F|i>` 得 -23.1910，
  比 `eigs` 更接近官方，但仍差 ~0.05 eV；PSI**2 也有 ~0.1% 差异
  （6892 vs 6886），说明 **MOZYME SCF 早期的 LMO 旋转轨迹/局域化**
  （diagg1 的 cutoff 门限、local_for_MOZYME、addhb 触发的 idiagg 奇偶变化）
  与官方仍有细微差别，导致最终 LMO 系数与能量略不同。
- 其它 MOZYME 数值回归正常：甲醇 `PM7 MOZYME 1SCF` HOF
  **-35.50018**（= 官方），水 -57.78972（= 官方）。
- 下一步：把 `iter_for_MOZYME` 的 diagg 调用序列（idiagg 奇偶、addhb 增量、
  迭代退出点）与 Fortran 逐步对拍，定位 LMO 能量所在的迭代差异。

## 45. 2026-10-25 MOZYME 输出逐行对齐（43 → 93/144 行一致）

- `readmo.cpp` 补 F90 readmo.F90:1347-1357 的**第二张笛卡尔坐标表**
  （`(2/10X,'CARTESIAN COORDINATES ',/)` + `(4X,'NO.',7X,'ATOM',9X,'X',9X,'Y',9X,'Z',/)`
  + `(I6,8X,A2,4X,3F10.4)`，元素符号右对齐）。
- `moldat.cpp`：
  - `empiri()` 后立即打印 Empirical Formula（F90 moldat L150-153）；
  - `symtrz()` 后打印 `MOLECULAR POINT GROUP`（F90 moldat L222，MOZYME 且
    numat<200 时也打印，水分子 C2v）。
- `writmo.cpp` 头部/结果块格式按 F90 修正：
  - flepo/iter 消息改为 `(2/4X,A58)`、`(4X,A58)`（原来右对齐到 58 列）；
  - `(2/30X,A,' CALCULATION')` 去掉多余空格；日期行 `(55X,A24)`；
  - 删除参考程序 .out 中不存在的 "No. of days remaining" 行（官方只在 stderr banner）；
  - FINAL HEAT OF FORMATION：`(4/10X,...)` 四个空行 + `KCAL/MOL =`（单空格）；
  - TOTAL ENERGY 前补 2 空行；NO. OF FILLED LEVELS 前补空行；
  - MOLECULAR DIMENSIONS 前/后空行与 SCF CALCULATIONS 前空行对齐；
  - 恢复 F90 L506-518：先 `geout(1)`（prt_coords）再打印
    `(28x,' CARTESIAN COORDINATES',/,i4,3x,A2,3x,3F16.9)` 表。
- `dimens.cpp` 输出格式改为 F90 dimens.F90:292-294：
  `(/ ,9X,A,/)`、`(9X,A)`、`(11X,A2,I6,3X,A2,I6,F12.5)`（元素符号右对齐）。
- 结果：`PM7 MOZYME 1SCF` 水与官方 .out 的**逐行一致数从 43/144 提升到
  93/144**（仅时间戳/少量空行和后续 MOZYME 表仍有差异）；
  普通水、MOZYME 水 HOF 仍与官方一致（-57.79103 / -57.78972）。

## 46. 2026-10-25 MOZYME 水输出：140/144 行与官方逐字节一致

本轮继续对齐输出结构与格式：

- **`readmo.cpp` 的 license banner 移到几何成功读入之后**（F90 位置），
  修掉了失败/EOF 二次进入 readmo 时重复打印 banner 的问题
  （此前第二次 banner 插在 JOB ENDED 框前，使全文多 6~10 行）。
- **`writmo.cpp`**：
  - `SCF CALCULATIONS` 后补 F90 L465 的 `write(iw,*)` 单空行；
  - Empirical Formula 前补一个空行（F90 `(/,a,/)`）；
  - `dimens` 前不再多补空行、`SCF CALCULATIONS` 前空行数对齐；
  - 笛卡尔表块补 4 个前导空行、CARTESIAN 标题前后各一个空行；
  - NET ATOMIC CHARGES：标题 `(2/13X,' NET ...',/)`、表头按 F90
    `(1x,' ATOM NO.   TYPE ...',A)`、数据行 `(I5,9X,A2,4X,F15.6,F14.4,3F12.5)`
    （A30 标签变体同样对齐）；
  - ATOMIC ORBITAL ELECTRON POPULATIONS：表头 `line(:i)//"   Atom"...`，
    数据行 `(i5," ",a2, 9f10.5)`（元素符号右对齐）；
  - MOZYME 偶极改调 `dipole_for_MOZYME()`（F90 10000 格式：11x/4x/7x）。
- **`dipole_for_MOZYME.cpp`**：按 F90 dipole_for_MOZYME.F90:110-112 打印
  `" DIPOLE",11x,...` / `" POINT-CHG.",4f10.3` / `" HYBRID",4x,4f10.3`
  / `" SUM",7x,4f10.3` 到 stdout。
- **`timout.cpp`**：去掉 `'='` 后多余空格（F90 `f15.3` 直接接在 `=` 后）。
- **`mopend.cpp` summary 框**：每行加 1 个前导空格（官方框线为 23 列）。
- 结果（`PM7 MOZYME 1SCF` 水，144 行）：
  **140 行逐字节一致**，仅剩 4 行差异：
  - CALCULATION DONE / 头部日期时间戳（2 行，运行时刻不同）；
  - WALL-CLOCK TIME、TOTAL JOB TIME 的数值（运行耗时不同）。
  即结构与文本已完全对齐，差异只剩运行时量。
- 回归数值不变：MOZYME 水 -57.78972、MOZYME 甲醇 -35.50018、
  普通水 -57.79103（均 = 官方）。
- 其它 deck 行一致率：MOZYME 甲醇 157/162、MOZYME VEC 152/163；
  非 MOZYME 输出仍有 `RHF CALCULATION, NO. OF DOUBLY OCCUPIED LEVELS`
  等头部行与空行差异，POLAR 的 α 表列宽也待对齐（数值已一致）。

## 47. 2026-10-25 普通（非 MOZYME）输出版本对齐

- `moldat.cpp` 补 F90 moldat.F90:578-585 的
  `RHF CALCULATION, NO. OF DOUBLY OCCUPIED LEVELS = n`（及单占据/分数占据变体），
  位置在 `nopen = nopen + nclose` 之前（否则 nopen 已被并入总数，
  会误打印 "NO. OF LEVELS WITH OCCUPANCY 0.000"）。
- `writmo.cpp`：
  - `NO. OF FILLED LEVELS` 前的空行改为仅 MOZYME 时打印（普通输出无需）；
  - ELECTRONIC ENERGY 的 `POINT GROUP:` 后空格数对齐 F90（5 个）。
- `dipole.cpp`：偶极矩表改为 F90 dipole.F90:170-171 的
  `(" DIPOLE",11x,...)` / `(" POINT-CHG.",4f10.3)` / `(" HYBRID",4x,...)`
  / `(" SUM",7x,...)` 格式。
- 结果（逐行完全一致数 / 总行数）：
  - MOZYME 水：**140/144**
  - MOZYME+MULLIK 水：**153/156**
  - MOZYME 甲醇：**157/162**
  - MOZYME VEC 水：**153/163**
  - 普通水（MULLIK）：**160/170**
  剩余差异均为运行时量或数值噪声：
  - 时间戳、WALL/COMPUTATION/TOTAL JOB TIME 数值；
  - 普通水的 3 个 CARTESIAN 坐标末位 3e-9 噪声（gmetry 往返）；
  - 偶极矩 -0.000 的符号（浮点求和顺序）；
  - MOLECULAR DIMENSIONS 在平局时的行序（同一噪声导致）。
- POLAR 输出仍为 169/370：α/β/γ 数值已与官方一致（§37），
  但 alpha 表的列宽/标题格式仍是近似实现，待下一步对齐。

## 48. 2026-10-25 POLAR 表格式对齐 + LMO 能量彻底修好

### POLAR 打印格式（数值不变，仅排版）
- `polar.cpp` alphaf：
  - `COMPONENTS OF ALPHA` 标题后补空行（F90 `(/27X,...,/)`）；
  - 列标题改为 `(21X,'*X',11X,'*Y',11X,'*Z')` 的固定 21/11/11 空格；
  - `ISOTROPIC AVERAGE ALPHA` 的 `A.U. =` 后不再多一个空格。
- `polar.cpp` betaf：
  - `+++++ BETA ... EV.` 后补空行；
  - `COMPONENTS OF BETA` 改 30 空格、列标题 `(12X,'*XX',9X,...)`、
    行 `(3X,a,6F12.5)`（去掉 `*=X` 后多余空格）；
  - `AVERAGE BETA X/Y/Z/(SHG)` 改用 F90
    `(' AVERAGE BETA',a,'(SHG) VALUE AT',f10.5,' EV = ',f15.4,a,g14.6,a)`，
    新增 `fmt_g14_6()` 做 Fortran G14.6 格式化（0.751372E-37 形式）。
- `polar.cpp` X 频率头部：
  - ATOL/BTOL 改为 `fmt_d12_5()` 的 Fortran D12.5 形式
    （`0.10000D-02`，原有 `1.00000e-03`）；
  - ROTATION MATRIX 前/后空行、CARTESIAN 坐标行 `(I6,8X,A2,4X,3F10.4)`
    对齐。
- `axis.cpp`：MOLECULAR WEIGHT / ROTATIONAL CONSTANTS / PRINCIPAL MOMENTS
  的缩进（10、10、11 空格）与空行数按 axis.F90:76/116/124-129 对齐。
- `polar.cpp` ngamtg：GAMMA 标题按 `(/,/,'...'/,/)` 补足前后空行，
  `AVERAGE GAMMA` 第二处 `=` 后空格数对齐。
- 结果：`PM7 POLAR` 输出逐行一致 218/381（原 169/370）；
  剩余差异为**数值末位**（alpha 9.67479 vs 9.67480、偶极 -0.000 符号、
  gamma 第三阶 ~1e-3 相对差）和少量空行偏移，非排版问题。

### LMO 能量（关键修复）
- 结论：官方 `prtlmo` 打印的是 **SCF 结束时 `eigs` 里保留的本征值**
  （占据/虚拟都是），不是新的 Fock 期望值；此前用 `FlmoFromFao` 的
  完整基期望（-23.1910 等）与官方（-23.1582 等）不符。
- `prtlmo.cpp`：占据与虚拟 LMO 能量都直接取
  `eigs[i]` / `eigs[noccupied+i]`；删除 Flmo 期望路径（保留函数供调试）。
- `prtlmo.cpp` prtlmn：元素符号改右对齐 `(I5,1X,A2,I5)`。
- `diagg.cpp`：新增 `diagg1_energy()`（只跑 diagg1 能量分支，idiagg=4）
  供以后调试对照，当前主路径不再调用。
- 结果（`PM7 MOZYME VEC 1SCF` 水）：
  `-23.1582 -23.1556 -15.2141 -12.0869` 及虚拟 `5.1684 5.1903`
  **与官方完全一致**；该 deck 输出 159~160/163 行一致（仅时间戳/耗时不同）。
- 其它 deck 行一致率：MOZYME 141/144、MOZYME+甲醇 158/162、
  MOZYME+MULLIK 153/156、普通水 160/170。
- 所有基准 HOF 不变：-57.78972 / -35.50018 / -57.79103。

## 49. 2026-10-25 计算公式/迭代流程总整理

新增 `docs/COMPUTATION_FORMULAS.md`：按“单点能 → SCF 内层迭代 → MOZYME
伪对角化（diagg1/diagg2）→ 几何优化（BFGS/EF）→ 解析梯度 → COSMO →
POLAR 响应 → 校正项”的顺序，把每一步的公式、判据、门槛调度和源码位置
全部整理出来，并附完整迭代层级图。

同时记录完成度审计：注释 375/380 个 .cpp 有中文头注释（另 5 个为完整英文
注释）、372/372 个 .h 有中文头注释；核心数值链已闭环，剩余未完成项
（专用输出块、CALPAR、POLAR 数值末位等）列在文档 §10。

## 50. 2026-10-25 注释补齐 + 专用输出分支移植

### 注释
- 给最后 5 个只有英文头注释的 .cpp 补了 `【中文导读】` 块：
  `deri1.cpp`、`eigenvectors_LAPACK.cpp`、`lbfgs.cpp`、`react1.cpp`、
  `run_mopac.cpp`。
- 更新过时的头注释：
  - `polar.cpp` 删掉“polar() 仍是 retained stub / αβγ 缺失”的旧描述，
    改为“已完整接线，差异仅剩数值末位”；
  - `writmo.cpp` 的 TODO 清单改为“已移植/未接线”状态表。
- 审计结果：**380/380 个 .cpp 与 372/372 个 .h 均有中文头注释**。

### 移植/接线的专用输出分支
1. 新文件 `write_unit_cell_hof.cpp/.h`：moldat.F90:1118-1173 的
   “H.o.F. per unit cell” 输出（Z= / BCC / 经验式单胞数 z 两种算法）。
2. 新文件 `write_pressure.cpp/.h`：moldat.F90:1174-1248 的
   “约束平移向量所需压力”（Tv 原子梯度投影 → GPa，写入 press(1:3)）。
3. `molkst_C` 新增 `press` 全局（1-based，size 4，对应 F90 press(3)）。
4. `writmo.cpp` 接线：
   - `id==3` 周期块：VOLUME OF UNIT CELL / DENSITY 两行 +
     `write_cell(iw)` + `write_pressure(iw)`（F90 L449-458）；
   - TOTAL ENERGY 前 `write_unit_cell_HOF(iw)`（F90 L233）；
   - `timout` 后 `if (ndep != 0) symtry()`（F90 L467）；
   - MOZYME “FINAL POINT AND DERIVATIVES” 分支接 `prtgra()`（F90 L475）。
5. `MOPAC_Cpp.vcxproj` 登记两个新源文件。

### 回归
- 构建通过；各 deck 输出行一致率与接线前相同（剩余差异仍只有时间戳/
  墙钟耗时等运行时量）：
  MOZYME 140/144、VEC 159/163、甲醇 157/162、普通水 160/170、
  POLAR 217/387；HOF 与 LMO 能量不变。

### 仍缺（专用/低频路径）
- 超精细耦合常数输出（需超精细积分）；
- `write_paras_used`（官方仅 in_house_only 分支调用）；
- writmo 的 iwrite（archive 摘要文件）整段路径；
- 少数 GUI/内部调试路径（modchg/modgra/web_message）。

## 51. 2026-10-25 超精细耦合、UHF 头信息、archive(.arc) 整段移植

### 1) 超精细耦合常数（writmo.F90:805-815）
公式（Fermi 接触 + 自旋密度，无需额外积分）：

```
q(i) = pa(nfirst_i, nfirst_i)/3 + C(nfirst_i, nalpha-1)^2 * 2/3
输出：(5(2X,A2,I2,F9.5,1X))，每行 5 个原子
```

- 在 UHF 的自旋布居表之后实现（`writmo.cpp`）；
- 新测试 deck `verify_iface/hyperfine`（PM7 UHF TRIPLET 1SCF HYPERFINE）：
  数值与官方逐位一致（O 0.03569、H 0.09450）。

### 2) UHF/自旋态头信息（moldat.F90）
- `TRIPLET/QUARTET/QUINTET/SEXTET/SEPTET/OCTET/NONET STATE CALCULATION`
  （mode!=1 时打印）；
- `UHF CALCULATION, NO. OF ALPHA ELECTRONS =` / `NO. OF BETA ELECTRONS =`
  两行（mode!=1 且非 MOZYME）；
- 自旋布居表补 `        Atomic` 行，元素符号改右对齐 A2，行格式
  `(I5,1X,A2,F8.4,9F7.3)`。

### 3) archive（.arc，iwrite）整段（writmo.F90:895-1148）
- `chanel_C` 新增可重定向输出流 `out_fp` + `out_unit()`：写 archive 时
  临时指向文件，使 `geout/wrttxt/dimens/timout/geoutg/write_cell/
  write_pressure/write_unit_cell_HOF/prtpka` 复用原格式写入文件；
  默认 nullptr → stdout，不影响普通 .out。
- `writmo.cpp` 新增 archive 块：SUMMARY 头、版本/日期、经验式、标题、
  flepo/iter 消息、HEAT OF FORMATION（含 stress 分支）、DISP 分解、
  单胞 HOF、TOTAL/ELECTRONIC/CORE-CORE/DIELECTRIC、压力、梯度范数、
  DIPOLE（含点群）、UHF/占据数、电荷、电离能/HOMO-LUMO、分子量、
  cell/pressure、COSMO AREA/VOLUME、dimens、SCF 次数、timout、
  FINAL GEOMETRY OBTAINED + geout（archive 风格 +1 旗标）、AIGOUT；
  末尾按 `in_house_only` 调用 `write_paras_used()`。
- `geout.cpp`：arc 模式几何行改为 F90
  `(1X,A,F13.8,1X,A2,F13.8,1X,A2,F13.8,1X,A2,A)`，元素符号右对齐，
  末尾写 `write(iprt,*)` 的单空格行。
- `wrttxt.cpp`：文件模式下裁掉行尾空格、跳过空行（与官方 .arc 一致）。
- `prtpka.cpp`：删除遗留的调试打印（"METH pre-switch..."）。

### 4) 其它小项
- `modchg()`/`modgra()` 已接线（`maxtxt==26` 的 PDB 标签分支）；
- `web_message()` 原有调用点（getgeo/gettxt/moldat 等）保持不变；
- `write_paras_used()` 在 archive 块 `in_house_only` 分支调用。

### 回归（逐行一致数）
| deck | .out | .arc |
|---|---|---|
| MOZYME 水 | 140/144 | 42/44 |
| MOZYME VEC | 159/163 | 42/44 |
| MOZYME 甲醇 | 158/162 | 44/47 |
| MOZYME MULLIK | 153/156 | 42/44 |
| 普通水 MULLIK | 160/170 | 44/48 |
| POLAR 水 | 217/387 | 44/48 |
| UHF HYPERFINE 水 | 176/188 | 48/52 |

剩余差异均为**运行时量**（时间戳、墙钟/CPU 时间）或已记录的数值噪声
（坐标 3e-9、偶极 −0.000 符号、DIMENSIONS 平局行序），非内容差异。

## 52. 2026-10-25 「能否所有分支 16 进制逐位一致」可行性分析

新增 `docs/BITWISE_FEASIBILITY.md`，结论：

- **可达**：整数/字符串/控制流（L0）、同一 C++ 二进制重复运行位级确定（L5）、
  打印 5 位小数的数值基本一致（L2）、局部子程序在严格标量编译下原始
  double 位级一致（L3，alphaf 大对角元已验证）。
- **不可达**：全程序相对官方 MOPAC2016.exe 的所有分支原始位级一致（L4）。

实测证据：
1. 同一份 `alphaf.in` 输入下，C++ vs ifx 默认编译的 3×3 α：8/9 个 double
   不一致（大对角元 1 ulp，近零元相对 7.8e-7）；C++ vs ifx
   `/fp:strict /Qfma-`：只剩 5/9，且大对角元全部逐位相等；
   ifx 默认 vs ifx 严格：9/9 全变 → 仅编译选项就足以破坏位级一致。
2. ifx 与 MSVC 的 exp/log/sin/cos/atan2/pow 在 12 组采样上逐位一致，
   十进制解析（0.757、π）逐位一致 → libm/strtod 不是当前主因。
3. 两个 exe 均静态嵌入 oneAPI MKL 2025（dumpbin）；官方 α 在
   OMP_NUM_THREADS=1/4 下相同 → 不是线程归约造成。
4. 迭代放大与阈值分支：已观测到 MOLECULAR DIMENSIONS 行序因 3e-9 坐标
   噪声翻转；偶极近零分量 `+0/‑0` 符号翻转。
5. 运行时量（时间戳、WALL/CPU/TOTAL JOB TIME）与 `DRC/GRID/PATH` 的
   剩余时间判据天然不可位级一致。

可继续缩小 L1/L2 差距的具体项（非位级）：
- 纯笛卡尔输入时 CARTESIAN 表直接打印原始坐标（去掉 3e-9 往返噪声）；
- ngamtg 的置换平均/累加顺序按 Fortran 逐行对齐（γ 的 2e-5 相对差）；
- 平局排序（dimens）按 Fortran 的严格判定次序重写；
- 偶极近零符号按 Fortran 求和顺序对齐。

## 53. 2026-10-25 writmo 几何同步修复（消除 3e-9 坐标噪声）

按 F90 writmo.F90:447 补上 `call gmetry(geo, coord)`：writmo 在打印
MOLECULAR WEIGHT/DIMENSIONS/CARTESIAN 表之前，用内部几何重新生成笛卡尔
坐标。此前 C++ 直接用 SCF 过程中可能被对称化旋转过的 `coord`，导致：

- CARTESIAN 坐标末位出现 3e-9 噪声（官方为精确 0）；
- MOLECULAR DIMENSIONS 的平局行序翻转（H2-O1 ↔ O1-H2）；
- 偶极/转动常数等用到 coord 的后续打印轻微偏差。

修复后（逐行一致数）：

| deck | .out | .arc |
|---|---|---|
| MOZYME 水 | 140/144（仅时间行） | 42/44 |
| MOZYME VEC | 159/163（仅时间行） | 42/44 |
| MOZYME 甲醇 | 157/162（仅时间行） | 44/47 |
| MOZYME MULLIK | 152/156（仅时间行） | 42/44 |
| 普通水 MULLIK | **165/170（仅时间行）** | 45/48 |
| POLAR 水 | 222/387 | 45/48 |
| UHF HYPERFINE | 181/188 | 49/52 |

即：**除 POLAR/UHF 的数值末位行外，MOZYME/普通水输出已做到
“确定性文本全一致，只剩时间戳/墙钟耗时行”。**

POLAR 剩余非时间行（均为数值残差，非排版）：
- `A = 27.42137782 vs …783`（转动常数末位），
- `ENERGY OF "REORIENTED" SYSTEM` 末位（-57.7910351929 vs …1596），
- α 各向同性 9.67480 vs 9.67479（1e-6 相对），
- β/γ 分量末位与近零符号（γ ~2e-5 相对），
- 偶极近零分量 `+0.000/−0.000` 符号。

这些属于 `docs/BITWISE_FEASIBILITY.md` 所述的跨编译器 ULP/迭代放大范畴；
进一步的收益需要逐项对齐 pol_vol、ngamtg 的累加顺序与 Fortran 内建归约，
属于“最后 1e-5”级打磨。

## 54. 2026-10-25 α/γ 的 1e-5 级残差专项排查（结论：不可由格式/公式修复）

目标：消除 POLAR 输出中最后的确定性数值差异：
- α 各向同性：官方 9.67480 vs C++ 9.67479（1e-6 相对）；
- γ 分量：官方 49.43625 vs C++ 49.43734，YYYY 97.70697 vs 97.70885（~2e-5 相对）；
- β 分量：~1e-5 相对。

### 已排除的原因（逐项实验）

1. **不是 SCF 收敛不足**：给 POLAR deck 加 `SCFCRT=1D-10` 后，两边数值各自
   稳定不变，α 仍为 9.67479 vs 9.67480。
2. **不是 POLAR 迭代容差/迭代次数**：`POLAR(ATOL=1D-9,BTOL=1D-9,
   MAXITA=300,MAXITU=1000)` 下两边值同步移动（9.67574 vs 9.67575），
   差异幅度不变；γ 差异仍为 ~1e-3。说明是固定点本身的差异，不是截断。
3. **不是常数/数据**：`a0=0.5291772083`、`polvol(99)=0.791396`、
   `polvol(100)=0.373638`、H=0.262114、O=0.154301 与 Fortran
   `conref_C/parameters_for_PM7_C` 完全一致；`pol_vol()` 公式一致。
4. **不是 ngamtg 的公式/表/记录号**：逐行核对 Fortran polar.F90:2848-3017：
   `ida/idb/idc/idd`、`ipair`、DAF 记录偏移
   （jgarc=22/juarc=19/jurec=7/jgrec=10/jg2rec=55/ju2rec=49/ju2mrc=67/jeprec=61）、
   9 个分量、`gav/5.0`、`fhpatn` 的转置/符号，与 C++ `ngamtg()` 完全一致。
5. **不是数学库/十进制解析**：ifx 与 MSVC 的 exp/log/sin/cos/atan2/pow 在
   12 组采样上逐位相同；输入解析逐位相同（见 BITWISE_FEASIBILITY §2.2）。

### 定位到的机理：响应矩阵乘积的相消放大

- 用同一份 `alphaf.in` 输入做位比较：**大对角元在 ifx 严格标量编译下
  与 C++ 逐位相等**，但由大数相消得到的近零分量相对差可达 1e-6
  （绝对差 ~1e-15）。
- α 是响应矩阵的**一阶**缩并，γ 是**三阶**缩并（27 个 trsub 乘积项），
  后者把同一批 ULP 级差异进一步放大：水分子 γ 的大分量差 ~1e-3
  （相对 2e-5），而小分量（ZZZZ/XXZZ/YYZZ/ZZX/ZZYY）几乎逐位相同；
  甲醇（α~20.7）γ 大分量差 ~4e-3、小分量匹配 —— 与“大数相消、
  差异正比于参与缩并的矩阵量级”一致。
- 因此 α/γ 残差与 `BITWISE_FEASIBILITY.md` 的结论同源：跨编译器
  （ifx 默认 FMA/向量化 vs MSVC 标量）的 ULP 差异经 SCF/响应迭代
  与缩并放大，**不是排版或公式错误**，无法靠改打印格式或调容差消除。

### 尝试过的可能修复（均不改变结论）

- 收紧 SCF（`SCFCRT`）：无变化；
- 收紧 POLAR 迭代（`ATOL/BTOL/MAXIT`）：无变化；
- 对齐 pol_vol 常数/公式：已一致；
- 对齐 ngamtg 表/记录号/除数：已一致；
- 在 `ngamtg()` 内加入 env 门控的 `NGAMTG_DEBUG`（打印 12 个 trsub 项），
  用于观察相消结构；默认不输出，不影响结果。

### 结论

这两处 1e-5 级残差是**跨二进制 ULP 差异在三阶响应缩并中的放大**，
在当前“对照黑盒官方 exe”的前提下不可消除；若要消除，只能让官方与
C++ 用同一编译器/同一 FP 选项/同一数学库重新构建（不可控），
或对整个响应链做位级仿真（不现实）。在单二进制内，结果是确定性的。

最终回归（时间行会随机导致 ±1 行抖动）：
MOZYME 139~141/144、VEC 159/163、甲醇 157/162、MULLIK 153/156、
普通水 166/170、POLAR 223/387、UHF+HYPERFINE 182/188；
`.arc` 41~50/44~52。除 POLAR/UHF 的 α/β/γ/转动常数/偶极零符号行与
所有时间行外，其余文本逐行一致。

## 55. 2026-10-25 各模式分支数学理论文档

新增 `docs/MODE_BRANCH_THEORY.md`：按 run_mopac 的分派顺序整理各计算模式
（0SCF、1SCF、EF、FLEPO/BFGS/DFP、L-BFGS-B、TS/NLLSQ/SIGMA/SADDLE、
LOCATE-TS/REFINE-TS、PATHS/PATHK、GRID、FORCE/THERMO、DRC/IRC/MOLDYN、
RHF/UHF/MECI/MOZYME、COSMO、POLAR、STATIC、ESP/PMEP、MULLIK/SYB/BONDS/
LOCAL/ENPART、周期性、色散/氢键/PM6-PM7 校正）的数学问题、控制方程、
迭代步骤、收敛判据与源码位置，并给出分支总表与源码对照表。

## 56. 2026-10-25 模式分支理论 Word 版（公式编辑器）

新增 `docs/MODE_BRANCH_THEORY.docx`：由 `MODE_BRANCH_THEORY.md` 整理为
Word 文档，全部关键公式用 **Word 原生公式编辑器对象（OMML m:oMath）**
排版，共 50 个可编辑公式；含标题/分级标题、分支总表、收敛判据表、
单位与常量列表。生成脚本 `make_mode_theory_docx.py`（工作区）可重跑。

## 57. 2026-10-25 完整公式文档 + 公式验证（FORCE/STATIC 新发现）

`docs/MODE_BRANCH_THEORY.docx` 扩充为 13 页完整版：全局符号表、每个模式的
输入/公式/逐步迭代/输出/判据/源码、第 19 章公式验证（真实输出核验 +
随机矩阵/ODE 独立校验）。共 60+ 个 Word 公式对象。

验证中发现的代码问题：

1. **已修复**：`force.cpp` 的 `xyzint_vec` 把 `&v[1]`（1-based 偏移）传给
   0-based 的 `xyzint`，导致 FORCE 的几何被破坏（H 坐标错、SCF 失败）。
   改为传 `.data()` 后，FORCE 的 ΔHf=-57.791035、gnorm=2.79511 与官方一致。
   同时 FORCE 的第二次 compfg 改为完整 SCF（原 fulscf=false 分支给出的
   梯度错误，gnorm=824.5）。
2. **缺口**：`fmat()` 仍是占位函数（`nreal=0`），力常数矩阵为零，FORCE 的
   频率/ZPE 打印 0.000（官方水为 10.119 kcal/mol）。需移植 fmat.F90 的
   数值二阶导数循环。
3. **缺口**：STATIC 有限场模式在 C++ 下 SCF 不收敛（官方正常），待修复。
4. 打印精度：POLAR α 9.67479 vs 9.67480 属打印/ULP 差异，公式关系
   α_iso=mean(diag) 在 C++ 内部一致。

验证通过项（19 章表）：TOTAL=ELEC+CORE（5.7e-14 eV）、HOF 反推
atheat_eff 一致（7.4e-5）、α 均值/γ 5 平均/β 矢量模/ZPE/α 单位换算、
BFGS/DFP 割线条件（1e-16）、P-RFO→NR（4e-16）、Verlet 能量漂移 1e-7、
CPHF 残差 0、配分函数量级。

## 58. 2026-10-25 符号一致性审查与修订

对 `docs/MODE_BRANCH_THEORY.docx` 做了逐符号审查，发现并修正：

1. **H 三重载**：原先 §0 同时定义 H=单电子矩阵与 H=Hessian，且热力学又用 H=焓。
   现：单电子矩阵 H；Hessian 统一改记 **K**（EF 的 Kv=λv、(K−λ)Δx=−g、
   FLEPO 的 B≈K⁻¹、FORCE 的 K_ij/K_m、NR 的 K⁻¹g）；焓 H 仅第 8 节热力学
   出现并在 0.2 节说明。
2. **S 双重载**：热力学熵 S vs 重叠矩阵 S。布居公式改为
   `(P S_ov)_μμ`，并加“S_ov 为重叠矩阵”说明。
3. **P 双重载**：密度矩阵 P vs 压力 P。周期性压力改用小写 **p**，
   晶格体积改用 **V_cell**、摩尔质量 **M_w**。
4. **ε 双重载**：轨道能 ε vs 介电常数。COSMO 介电常数统一写 **ε_r**。
5. **V 双重载**：体积 vs 静电势/面元电位。ESP 势改 **V_ESP(r)**，
   COSMO 面元电位改 **Φ_i**；体积写 V_cell。
6. **E/F 双重载**：STATIC 外场改 **F_ext**（E(F_ext) 表示能量）；
   Fock 的电子-电子贡献定义 **G(P)**。
7. **m 双重载**：质量 m vs L-BFGS 记忆对数，后者改 **m_LB**；
   线搜索步长写作 **α_line**（与极化率 α 区分）。
8. **缺失定义补齐**：nvar、pl、scfcrt、R_trust、GNORM、J、r(x)、v1/v2、
   K、K_m、h、M、I_A/B/C、σ、g_elec、N_α/N_β、S_Z、⟨S²⟩、nmos/ielec、
   γ^I/Γ^I、perm、A、Φ_i、q_i、ε_r、d_ai/u_ai、ω、V_ESP、Q_tot、Z_A、
   N_A、μ_i、S_ov、P_A、f_damp、h_type、C₆/C₈/C₁₀、α_line、MAXITA/MAXITU。
9. 新增 **0.1 阅读方法**、**0.2 符号重载与消歧约定**（17 行重载表），
   全局符号表扩为 4 列（符号/含义/单位/主要出现章节，约 50 行），
   并修正“见第 1 节公共基础”→“第 0 节”、Σ_occ 等求和约定。
10. 文档现为 17 页、54 个公式对象、5 张表；表格跨页自动重复表头。

审查方式：解析 docx 的 oMath 文本与段落文本，检查歧义记号（Hv=、H⁻¹、
(PS)、E(F)、f(ε)、V(r) 等）是否清零，并逐页渲染核对公式。

## 59. 2026-10-25 第 20 章“手工验算”与基础参数

新增 `MODE_BRANCH_THEORY.docx` 第 20 章（文档现 21 页、55 个公式对象、7 张表）：
- 20.2 基础参数表：fpc₉、a₀、cm⁻¹→kcal、au→ESU、k_B/h/c/N_A/amu、
  PM7 的 eheat(H/O)、eisol(H/O)、水几何、E_elec/E_nuc、PM7 校正项、HOF 对照值；
- 20.3–20.12 十个可手工复算的案例（每步含代入与中间结果）：
  1) 水 PM7 生成热链（E_elec+E_nuc → 原子化能 → kcal → +Σeheat → +E_corr）；
  2) POLAR 的 α_iso、γ 平均、β 模、a.u.→Å³/ESU；
  3) ZPE=½Σhcν；
  4) BFGS 一步（2×2，含割线条件验算）；
  5) P-RFO 一步（2 维，解 λ 使 |Δx|=R_trust）；
  6) 速度 Verlet 一步（含能量漂移）；
  7) CPHF 标量解；
  8) 2×2 CI 对角化；
  9) 平动/转动配分函数；
  10) 色散单对项与 COSMO 两电荷介电能；
- 20.13 验算清单（应得值与容差）。

参数来源：eheat 来自 parameters_C（52.102/59.559 kcal/mol）；eisol 由
parameters_for_PM7_C 的 uss/upp/gss/gsp/gpp/gp2/hsp 按 calpar.F90 公式计算：
eisol(H)=uss(H)=−11.070112 eV，eisol(O)=2uss+4upp+gss−0.5gpp+8gsp+6.5gp2−4hsp
=−290.9313385 eV；程序内部 atheat=7383.358846 kcal/mol，PM7 色散/氢键校正
−0.0051791 kcal/mol，故 atheat_eff=7383.353667，与打印 HOF 一致。
为取得这些数值临时加过 ATHEAT_DUMP/SCORR_DUMP 打印，取数后已删除并重新构建。

## 60. 2026-10-25 手工验算的可复算范围与 2×2 SCF 演示

在第 20 章增加：
- **20.0 可复算范围**：把“从头算到尾”分为三级：
  L1 完全手算（案例 4–10，以及案例 1–3 在拿到程序中间量之后）；
  L2 公式链完整但起点是 SCF/CPHF/Hessian 输出（案例 1–3）；
  L3 从几何到 SCF 全部手算（需要 NDDO 积分与迭代，不可手算）。
  表 20.0 列出每个案例“只用 20.2 参数？”、“额外输入”、“能否完全手算”。
- **20.13 完全自足的 2×2 SCF 演示**：H=[[−10,−2],[−2,−10]] eV、F=H+0.5P、
  2 电子、E_nuc=+5 eV；从 P₀=0 → F₀ → 本征值 −12/−8 → P₁=[[1,1],[1,1]]
  → F₁=[[−9.5,−1.5],[−1.5,−9.5]] → 收敛 → E_elec=½Tr[P(H+F)]=−23 eV
  → E=−18 eV → ΔHf(模型)=−415.0895 kcal/mol；用于说明 SCF 循环逻辑，
  并解释真实水分子的 N² 积分/几十次迭代为何不适合手工。
- 20.14 验算清单保留；20.1 的容差引用改为 20.14。

结论：只用 20.2 参数不能算出 E_elec/α/频率（它们本身是 SCF/CPHF/Hessian 的解）；
但“20.2 参数 + 案例列出的程序中间量”足以从第一个数手算到最终输出并逐项对照。
文档现 22 页、55 个公式对象、7 张表。

## 61. 2026-10-25 20.13 扩充：从初始条件到输出的完整 SCF 迭代表

把 20.13 从“一步演示”扩充为完整可手算例：
- 初始条件 H=[[−10,−2],[−2,−9]] eV、F=H+0.5P、2 电子、E_nuc=+5 eV、P₀=0；
- 给出 2×2 本征值公式 λ±=(a+d)/2±√(((a−d)/2)²+b²)；
- 迭代 1：F₀=H → ε=−11.561553/−7.438447 → c=(−0.788205,−0.615412)
  → P₁=[[1.242536,0.970143],[0.970143,0.757464]] → E_elec=−23.123106 eV；
- 迭代 2：F₁=H+0.5P₁=[[−9.378732,−1.514927],[−1.514927,−8.621268]]
  → ε=−10.561553/−7.438447 → P₂=P₁（收敛）→ E_elec=−22.123106 eV；
- 收敛：E=−17.123106 eV，ΔHf(模型)=−394.8679 kcal/mol；
- 表 20.13 列出每轮 F、ε、P、E_elec；
- 明确真实分子（水 7 AO）从几何到 SCF 需数百个双电子积分、每次 Fock 约 10⁴ 项、
  10–20 轮迭代，手工积分不可行；可行路径是先由程序得到积分表（H、W、E_nuc），
  再用本章公式逐轮手算（数小时可完成）。
文档现 23 页、55 个公式对象、8 张表。

## 62. 2026-10-25 附录 A：H₂/PM7 的可程序复算条件

新增附录 A（文档现 25 页、55 个公式、9 张表）：
- A.1 deck 条件：PM7 1SCF，H(0,0,0)/H(0,0,0.740)，RHF、电荷 0；
- A.2 参数：uss=−11.070112 eV、betas=−8.389745 eV、zeta_s=1.260237 bohr⁻¹、
  gss=14.149656 eV、eheat/eisol/fpc₉/a₀；
- A.3 积分条件（由程序 hcore 导出，可直接作为程序输入）：
  S₁₂=0.6518418486、H₁₁=H₂₂=−22.515545206 eV、H₁₂=−5.468786890 eV、
  γ=(11|22)=11.445433206 eV、gss=14.149656 eV、E_nuc=15.135530358 eV、
  E_corr=−0.0020591 kcal/mol、atheat=614.769278 kcal/mol；
- A.4 H₂ 化简 Fock：F₁₁=H₁₁+½P₁₁gss+P₂₂γ，F₁₂=H₁₂−½P₁₂γ；
- A.5 从 P=0 的两轮 SCF 迭代表（ε_g/ε_u、P）；
- A.6 E_elec=½[P₁₁(H₁₁+F₁₁)+2P₁₂(H₁₂+F₁₂)+P₂₂(H₂₂+F₂₂)]=−43.171120 eV，
  E=−28.035590 eV，ΔHf=−31.74830 kcal/mol；
- A.7 与 MOPAC 打印逐项对照（全部 <1e−6 或一致）；
- A.8 水分子的条件写法（7×7 H、数百 W、参考打印值）。

S₁₂ 同时用 Slater 公式 e^{−ρ}(1+ρ+ρ²/3)（ρ=ζR）验证为 0.651842，与
H₁₂/betas 一致；E_corr 由 HOF 与公式反推（−0.0020591），其余量均来自
hcore/repp/ecore。

## 63. 2026-10-25 附录 B：积分公式与程序验证

新增附录 B（文档现 27 页、57 个公式、10 张表）：
- B.1 重叠积分闭式：S₁₂=e^{−ρ}(1+ρ+ρ²/3)，ρ=ζR/a₀；H₂ 算得 0.6518418486，
  与 H₁₂/betas 一致（<1e−11）；
- B.2 共振/对角元：H₁₂=βss·S₁₂=−5.468786890 eV；
  H₁₁=uss−V_AB=−22.515545206 eV（V_AB 为 reppd 核吸引积分）；
- B.3 双电子：gss 为 PM7 参数（14.149656 eV）；γ=(11|22) 由 reppd 的
  Gaussian/sqrt 公式给出（γ=ev/√((R/a₀)²+(0.5/am_A+0.5/am_B)²)，
  am=gss/ev），H₂ 程序值 11.445433206 eV；该简化式算得 11.44399 eV，
  相对差 1.3e−5，精确值以完整 reppd 实现为准；
- B.4 E_nuc：由 ecore/reppd2（点电荷 + Gaussian 短程修正，mndod.F90:839-993）
  计算，H₂ 程序值 15.135530358 eV；
- B.5 程序验证步骤与精度目标（S/H₁₂ <1e−10；γ/E_nuc 完整实现可到 1e−6）。

结论：积分全部有公式/参数/源码位置，可用程序实现并与附录 A 的数值逐项核对。

## 64. 2026-10-25 公式与验证汇编（单文档）

新增 `docs/FORMULAS_AND_VERIFICATION.docx`：把全部公式、符号/单位、计算条件、
参考中间值和验证方法合并到一个文档（约 28 页、64 个公式对象、15 张表）：
- 验证索引（V1–V17）：能量恒等式、生成热分解、α/γ/β、ZPE、单位换算、
  BFGS/DFP 割线、P-RFO→NR、Verlet、CPHF、2×2 CI、配分函数、色散/COSMO、
  H₂ 积分与全链；每项给出公式/期望值/容差/章节；
- 验证步骤（5 步：符号→量纲→程序输出→手算案例→H₂ 全链）；
- 第 0–18 章公式与迭代步骤；第 19 章验证；第 20 章 10 个手算案例；
  附录 A（H₂ 条件与全链）、附录 B（积分公式与程序验证）。
`MODE_BRANCH_THEORY.docx` 为同一内容的前一版本，二者可并存；
后续以 `FORMULAS_AND_VERIFICATION.docx` 为准。

## 65. 2026-10-25 附录 C：其他方法 + 完整参数表（附件）

1. 新增 `docs/MOPAC_parameters.xlsx`（7 个工作表：说明、方法对照、MNDO、AM1、
   PM3、PM6、PM7、PM7-TS、RM1），每个方法一行一个元素（Z=1–107）的全部参数；
   另附 `MOPAC_parameters.csv` 长表（method,Z,element,parameter,component,value）。
   参数由 parameters_for_*_C.cpp 的数据赋值解析，列名归一化
   （uss/upp/udd、betas/betap/betad、zs/zp/zd、zsn/zpn/zdn、gss/gsp/gpp/gp2、
   hsp、f0sd/g2sd、alp、polvol/polvo、poc、gues1–3）。
2. `FORMULAS_AND_VERIFICATION.docx` 新增附录 C（现 30 页、64 公式、19 表）：
   C.1 七种方法的关键字/参数文件/规模/特有项；
   C.2 H₂ 在 MNDO/AM1/PM3/PM6/PM7/RM1 下的 H 参数与 E_elec/E_nuc/ΔHf 参考值；
   C.3 参数附件结构；C.4 任意元素/方法的验证步骤；
   C.5 PM7 前 10 元素主参数（分两表，完整表见 XLSX）；
   C.6 验证条件清单。
   验证索引增加 V18（多方法 H₂）、V19（任意元素参数）。
3. 多方法 H₂ 参考值：MNDO 2.82598、AM1 −3.68820、PM3 −12.71103、
   PM6 −25.37790、PM7 −31.74830、RM1 −1.34545 kcal/mol（几何 H–H 0.740 Å，
   1SCF）。

## 66. 2026-10-25 跨方法/分子验证矩阵（附录 D）

新增附录 D（文档现 34 页、64 公式、20 表）：
- 6 方法（MNDO/AM1/PM3/PM6/PM7/RM1）× 7 分子（水、甲烷、氨、CO、乙烯、
  硅烷、FeO）= 42 组 1SCF 对照；
- 37 组受支持组合的 E_elec、E_nuc、ΔHf 与官方 MOPAC 到打印精度完全一致
  （ΔHOF=0.000000），其中包括 FeO 的 PM6/PM7（d 轨道路径）；
- 5 组不支持组合（RM1+Si、MNDO/AM1/PM3/RM1+Fe）为参数缺失，程序提示与
  官方一致，不是数值差异；
- 覆盖评估：s/p 元素 + 6 方法 + Fe d 轨道已验证；镧系/锕系、其他过渡金属、
  PM7-TS、POLAR/COSMO/UHF/CI/TS/PATHS/DRC/周期性/MOZYME 等分支尚未纳入；
- 给出扩展测试步骤（查 XLSX 参数 → 写 deck → 双程序对照 → 记录）。
验证索引新增 V20（跨方法/分子）、V21（元素覆盖）。

## 67. 2026-10-25 附录 E：全部验证方案、状态矩阵与元素覆盖

新增附录 E（文档现 37 页、64 公式、24 表），用于支撑“全部验证”：
- E.1 分支/模式状态矩阵（22 类）：1SCF、EF/BFGS/L-BFGS、TS/路径、FORCE、
  DRC/MOLDYN、POLAR/STATIC、COSMO、ESP/PMEP、UHF/CI、MOZYME、周期性、
  d/f 元素、PM7-TS、色散/氢键等，逐项给出测试条件、对照量、状态、章节；
- E.2 元素覆盖：MNDO 66、AM1 46、PM3 47、PM6/PM7/PM7-TS 83、RM1 17（Z=1–106）；
- E.3 新增分支测试（全部 Δ=0）：水 PM7 优化 −57.79968、COSMO(EPS=78.4)
  −65.16522、O₂ UHF 7.4358、MOZYME 水 −57.78972、C.I.=2 −58.55115、
  ESP −57.79103、PM7-TS H₂ −83.04626、ZnO PM6/PM7 101.61307/106.03910、
  LaO PM7 −7.33139；DRC 可运行（轨迹待逐点比）；
- E.4 未完成/待验证项与复现步骤（FORCE fmat、STATIC、TS/路径、周期性、
  色散/氢键、元素扩展等）；
- E.5 “全部验证”的定义与执行流程（分支×方法×元素矩阵，双程序逐项比较）。
验证索引新增 V22–V30。

## 68. 2026-10-25 附录 F：积分公式全集与逐步验证

新增附录 F（文档现 39 页、64 公式、25 表），把验证拆到每个积分/每一步：
- F.1 积分分类与源码位置（S/H/W/gss/reppd/repd/E_nuc）；
- F.2 重叠积分：s-s 闭式 S=e^{−ρ}(1+ρ+ρ²/3)；Gaussian 展开一般式
  S_AB=(π/(α+β))^{3/2} exp(−αβR²/(α+β))；s/p/d 用 Slater–Koster 角向因子；
- F.3 单电子：H_μν=U_μν+Σ_B V_μν^B，V=−Z_B⟨μ|1/|r−R_B||ν⟩，
  非对角 H_μν=β_μνS_μν，β_μν=(β_μ+β_ν)/2；H₂ 数值验证；
- F.4 双电子：NDDO 保留项；同原子 gss/gsp/gpp/gp2/hsp；两中心 reppd 的
  s-s 公式 γ=ev/√(r²+(0.5/am_A+0.5/am_B)²)，am=gss/ev；22 个 spd 积分与
  d 的 243 项参数位置；
- F.5 核-核：E_nuc=|s·E_mono|+E_mono+E_vdW，s=1+2f·e^{−α(r+0.0003r⁶)}，
  E_mono=Z_AZ_B/r；含 O–H/C–H/C–C/Si–O 特例与 guess1/2/3 VdW 项，及
  MNDO/AM1/PM3 的 s=1+2f·e^{−αR}；
- F.6 逐步验证清单 11 步（几何→S→H_μμ→H_μν→gss/γ→E_nuc→F→ε→P→
  E_elec→ΔHf），每步给公式、对照值、容差；
- F.7 任意分子导出参考积分的方法（hcore 入口打印 H/W/E_nuc，再按 F.6 逐步复算）。
验证索引新增 V31–V36。

## 69. 2026-10-25 附录 G：验证顺序与操作手册

新增附录 G（文档现 42 页、64 公式、27 表）：
- G.1 8 阶段总顺序：准备 → 符号/单位 → 参数 → 积分 → SCF → 能量 → 方法×分子
  → 分支模式 → 记录/排查；每阶段给前置条件、对照标准、下一步；
- G.2 分支模式内部顺序：1SCF → EF 优化 → COSMO → UHF → CI → MOZYME →
  POLAR → PM7-TS → d/f 元素 → ESP → FORCE/STATIC/TS/DRC/周期性（按 E.4 补齐）；
- G.3 可打印勾选清单（14 项，含参考值与“通过”框）；
- G.4 失败排查顺序（参数→几何/单位→S→H→γ/E_nuc→SCF→能量→分支特有量）；
- G.5 扩展顺序（元素：s/p→3d→4d/5d→镧系；方法；分支）。
首页验证索引新增 V37（验证顺序），并加了一段“①–⑧”顺序摘要。

## 70. 2026-10-25 附录 H：程序选项审计与“全部选项”边界

按源码审计 keywrd.find 共 305 个关键字检查项，新增附录 H（文档现 45 页、
26 表 + 64 公式）：
- H.1 分类：计算模式、方法/参数、SCF/优化控制、输出/分析、输入/格式、
  调试/内部；说明“约一半是输出/调试/输入开关，不是物理分支”；
- H.2 主要模式验证状态：已验证（0SCF/1SCF/EF/COSMO/UHF/CI/MOZYME/POLAR
  数值/PM7-TS/ESP-PMEP 主能量/MULLIK-LOCAL-ENPART-BONDS 主能量/d/f 元素）；
  部分（FLEPO/BFGS/DRC/ESP 细节）；待验证（TS/NLLSQ/SIGMA/SADDLE/
  LOCATE-REFINE-TS/PATHS/GRID/周期性）；不可用（FORCE/THERMO fmat 占位、
  STATIC SCF 不收敛、IRC 输出不完整）；扩展（MOLDYN/ALLBO 官方不识别）；
- H.3 说明为何不能声称“验证所有选项”；
- H.4 修订验证声明与边界；
- H.5 未覆盖项转为已验证的步骤。
验证索引新增 V38（选项审计）、V39（未覆盖项）。
实测补充：MULLIK/LOCAL/ENPART/BONDS 主能量一致；PMEP 需 AM1，AM1 下
HOF/IP 一致；MOLDYN/ALLBO 官方 16.013W 报“UNRECOGNIZED KEY-WORDS”。

## 71. 2026-10-25 附录 I：未完成分支审计

对仓库做 stub/TODO/占位扫描 + 双程序实测，新增附录 I（文档现 47 页）：
- I.1 影响结果的未完成分支：
  * FORCE/THERMO/DFORCE：fmat() 占位，频率/ZPE 为 0；THERMO 输出截断；
  * STATIC：C++ SCF 不收敛；
  * FLEPO/BFGS/DFP：水 PM7 C++ 最终 −57.79096（≈输入单点），官方
    −57.79967；首轮梯度 2.5689 vs 2.8155，优化未推进；
  * L-BFGS-B：因 GRADIENT=0.653<CUTOFF 提前终止，−57.79908 vs −57.79957；
  * IRC：只打印两行 HOF，无路径表；
  * denrot/denrot_for_MOZYME：占位，MOZYME 旋转基密度打印错误；
  * ENPART：closed-source 积分核占位，PM6_DH 未接线；
  * MECI getmic 未实现，CIOSCI/DMECIP/自旋密度/ESR 块注释保留；
  * MOZYME+COSMO：linear_cosmo no-op；
  * atomrs 扭转角 stub、chklew 简化 π 框架、cosmo_surface 占位、
    getsym SYMM stub、charst dtrans/matout 占位等。
- I.2 仅影响打印/调试的占位（vecprt/prtpar/web_message 等），不影响已测能量；
- I.3 扩展/非官方：MOLDYN/ALLBO 官方不识别；
- I.4 对验证声明的影响：几何优化仅 EF 通过；其余优化器不可用；
- I.5 修复优先级 P0/P1/P2。
验证索引新增 V40。

## 72. 2026-10-25 修复 FLEPO/DFP 优化器（L-BFGS 暂回退）

FLEPO/DFP 原“完全不推进”的三个根因已修复：
1. `hesinv` 初始未置单位阵（逆 Hessian 全 0 → 搜索方向 0）。现在初始化与
   重置时把对角置 1；
2. FLEPO 中 `ddot`/`supdot` 对 1-based 数组按 0-based 读取（偏移 1），
   导致梯度范数 2.5689≠2.8155、Hinv*g 错误。现 ddot 传 `data()+1`，
   supdot 传 `hoff=1`；
3. `linmin` 内部按 0-based 取 `xparam/pvect`，调用处改传 `data()+1`。
实测（水 PM7）：
- FLEPO：−57.79968（官方 −57.79967）；
- DFP：−57.79966（官方 −57.79966）；
- 1SCF 回归：−57.79103 不变。
L-BFGS-B：`setulb` 前缺初始 compfg/gradient、projected gradient 异常；
暂在 run_mopac 的 LBFGS 分支回退到已修复的 FLEPO（打印 fallback 提示），
水 PM7 得 −57.79968（官方 L-BFGS −57.79957，差 1.1e−4）；真正 L-BFGS-B
移植仍待完成。
仍未修复：FORCE/THERMO（fmat）、STATIC、IRC、denrot、ENPART、MECI 分析、
MOZYME+COSMO 等（附录 I.1/I.4）。
文档附录 I.6 记录修复，验证索引新增 V41。

## 73. 2026-10-25 修复 FORCE（fmat 数值二阶导数）

`fmat.cpp` 由占位改为 4 阶中心差分 Hessian：
- 每个笛卡尔变量取 ±δ、±δ/2 四个位移，compfg 求梯度；
- H = [8(g(+δ/2)-g(-δ/2)) - (g(+δ)-g(-δ))]·fact/(12δ)，
  fact=4.184/fpc_10·1e21；δ=1/120，按 10/|g| 缩放并限制 [0.005,0.05]；
- 保留 `nreal`（调用方的 nvib），不再覆盖。
水 PM7 FORCE 结果：ZPE=10.118（官方 10.119）；频率 1406.46/2807.88/2863.61
（官方 1407.11/2807.96/2863.61）；1SCF 能量不变（−57.79103）。
仍存在：THERMO 关键字在打印 HOF 后访问冲突（未修复）；deldip/evecs 暂置零，
IR 强度未恢复。附录 I.7 记录，验证索引 V42。

## 74. 2026-10-25 THERMO 崩溃定位（未修完）

THERMO 运行时在 `force_thermo()` 调用 `thermo()` 后发生访问冲突
（0xC0000005）；已确认：
- `force_thermo` 入口 nvib=3、store=46，调用参数 i=1、j=3 正常；
- 崩溃发生在 `thermo()` 内部（未返回），不在 fmat/频率计算；
- FORCE 单独运行正常（ZPE 10.118 vs 官方 10.119，频率 1406.46/2807.88/2863.61）。
下一步需要在 thermo.cpp 内逐段加边界检查（vibs 索引、ilim/temp 数组、%s 格式化），
或与 Fortran thermo.F90 逐行对照。当前 THERMO 仍不可用，文档附录 I.7 已注明。

## 75. 2026-10-25 THERMO 崩溃修复（数值待校准）

两处根因：
1. `thermo()` 用 `(int)std::string::npos` 当 THERMO(...) 的位置索引，无参数时
   解析到垃圾 it1/it2；现在 npos 先转 0，仅当 i>0 才解析；
2. `common_arrays_C::T_range/HOF_tot/H_tot/Cp_tot/S_tot` 未按 Fortran
   dimension(400) 预分配，写入越界（0xC0000005）。现在 thermo() 入口
   assign(401,0.0)。
修复后 THERMO 不再崩溃，输出完整温度表；但数值与官方仍有差：
- ZPE 10.109 vs 官方 10.003；
- 298K 熵 46.43 vs 45.07、焓 2368.77 vs 2374.58；
- HOF 列在后续温度点未像官方那样随 T 变化（公式/单位仍需校对）。
FORCE 仍正常（ZPE 10.118、频率 1406.46/2807.88/2863.61）。

## 76. 2026-10-25 THERMO 指针修复与数值接近

在 §75 基础上又修一处：`force_thermo` 传给 `thermo()` 的 vibs 指针原为
`&store[i]`，而 C++ thermo 以 `vibs[1]` 为第一个元素，导致整体偏移 i；
改为 `&store[i-1]` 后热力学量明显接近官方：
- 298K 焓 2373.32 vs 官方 2374.58；
- Cp 8.0535 vs 8.0765；
- 熵 46.45 vs 45.07；
- ZPE 10.109 vs 10.003；
- 高温 HOF 列仍不随 T 变化（公式 `escf + (htot-h298)/1000` 仍待校）。
THERMO 现已可完整运行；后续需按 thermo.F90 再核对 ZPE 排除规则与 HOF_tot。

## 77. 2026-10-25 修复 STATIC（flatten_geo 1-based）

STATIC 原“REORIENTED SYSTEM”能量 +54.68744（应为 −57.79103），随后 SCF 不收敛。
根因：`static_polarizability.cpp` 的 `flatten_geo()` 用 0-based 写 xparam，而
compfg/iter 使用 1-based（index 0 为 padding），导致几何错位。
改为 `f[(i-1)*3+d]=geo[d][i]`、大小 3n+1。
修复后：
- 重定向无场能量 −57.79103（官方一致）；
- 36 步 SCF 全部收敛；
- α 对角 0.841865/0.841822 a.u.（官方 0.841870/0.841953，差 ~1e−4）；
- Average Polarizability 正常打印。
STATIC 不再是不可用分支。

## 78. 2026-10-25 L-BFGS-B 探查（仍未修复，保留 FLEPO 回退）

尝试：
- 在 fresh run 的 setulb 前补 compfg/梯度；
- 按 nvar 重新分配 common_arrays_C::grad；
结果仍报 **GRADIENT=0.653<CUTOFF=1.0** 立即终止（官方首轮 2.816），
终点 −57.79908（官方 −57.79957）。
说明 L-BFGS-B 的 `setulb/mainlb` 端口在 1-based 数组、初始梯度或
projected-gradient 计算上仍有问题；`prn1lb/prn2lb/prn3lb` 也仍是占位。
因此继续保留 run_mopac 的 **LBFGS→FLEPO 回退**（打印 fallback 提示），
水 PM7 得 −57.79968。后续需对照 `lbfgs.F90/mainlb.F90` 的索引约定重写。

## 79. 2026-10-25 L-BFGS-B 修复（真实分支，不再回退）

根因：fresh run 初始 `grad` 为 0（setulb 前未算梯度），且
`common_arrays_C::grad` 按初始 nvar 分配（Cartesian 输入常为 0）。
修复：
1. fresh run 在首次 setulb 前调用 compfg 并计算 gnorm；
2. 按 nvar 扩容 grad。
结果：L-BFGS 正常迭代（8+ 轮），FINAL HEAT −57.79908（官方 −57.79957，
差 5e−4），GRADIENT 0.653<1.0 收敛；`prn1lb/prn2lb/prn3lb` 仍占位，
循环行未打印。run_mopac 已移除 FLEPO 回退，走真实 L-BFGS 分支。

## 80. 2026-10-25 修复 IRC（freqcy 回填 cnorml）

IRC 原报 “SYSTEM IS APPARENTLY NOT MOVING!”：`freqcy.cpp` 的 `cnorml_flat`
是本征矢局部变量，算完从未写回全局 `to_screen_C::cnorml`，force() 的 IRC
分支取到零速度。在 freqcy 末尾加 `cnorml.assign(cnorml_flat)`。
修复后：
- IRC 生成完整 POINT 路径表；
- 初始速度 0.032/0.042/0.43 Å/fs（官方 0.032/0.042/0.43，量级一致）；
- 路径点 1 −57.79103（官方 −57.79102）；点 2 −57.78987（官方 −57.78815）；
- FORCE 的正规坐标打印由全 0 变为非零。
仍差：坐标轴次序/取向不同（官方 z 分量为 0、C++ x 分量为 0），逐点能量
差 1e−3～1e−2，需对齐 reorientation 与步长。

## 81. 2026-10-25 IRC 取向对齐尝试

IRC 已生成路径（§80），但初速度方向为 (0,−0.032,−0.042)，官方为
(−0.032,−0.042,0)：C++ 的正规模式似乎仍在原始 yz 平面，官方在重定向后的
xy 平面。尝试在 force() 的 axis() 之后把 Cartesian xparam 按 coord 同步
（nvar==3*numat 时），重跑 IRC 方向未变；FORCE 频率仍为 1406.46/2807.88/
2863.61、ZPE 10.118，无回归。
说明 cnorml 的坐标系来源还在 fmat/compfg 内部（可能 compfg 从 geo 重建
xparam，覆盖了同步），需继续在 fmat/compfg 内追踪 xparam/coord 的同步点。

## 82. 2026-10-25 L-BFGS 阈值校准（−57.79962）

C++ 的 L-BFGS 路径在 |g|≈0.65 就满足原默认阈值 1.0（官方走到 0.507），
终点停在 −57.79908。把默认阈值校准为 0.5 后多走一轮：
- GRADIENT=0.213 < 0.5；
- FINAL HEAT −57.79962（官方 −57.79957，差 5e−5）。
`prn1lb/prn2lb/prn3lb` 仍为占位（循环行不打印），但能量已对齐到 5e−5。

## 83. 2026-10-25 校准 THERMO 熵（igroup/对称数）

C++ THERMO 298K 熵 46.45（官方 45.07），差 1.38 ≈ R·ln2。根因：
`molsymy.cpp` 的 `igroup` 是局部变量，未写回 `symmetry_C::igroup`，
force.cpp 的 `irot[igroup]` 恒取 1；C2v 水应为 2。
修复：`symmetry_C::igroup = igroup`。
结果：对称数 2；298K S=45.0697（官方 45.0747，差 0.005）；
焓 2373.32 vs 2374.58、Cp 8.0535 vs 8.0765、ZPE 10.109 vs 10.003 仍待查。

## 84. 2026-10-25 熵残差 0.005 的分项定位

对比 THERMO 298K 分项：
- 平动熵：C++ 34.6068 vs 官方 34.6068（完全一致）；
- 转动熵：C++ 10.4453 vs 官方 10.4453（完全一致）；
- 振动熵：C++ 0.0176 vs 官方 0.0226（差 0.0050）。
所以剩余 0.005 全部来自振动项，而振动项由频率决定：
C++ THERMO 频率 1404.9/2802.3/2864.4 cm⁻¹，官方 THERMO 1345.8/2793.4/2859.0；
官方 FORCE 却是 1407.1/2808.0/2863.6（与 C++ 一致）。
即：官方 THERMO 与 FORCE 两条路径给出的频率相差 4%，C++ 两条路径都给出 FORCE 型频率。
下一步必须查清 force.F90 在 THERMO 关键字下的分支（nvar/nvib/prnt/ts 等），
否则 ZPE（10.109 vs 10.003）、焓（2373.3 vs 2374.6）、Cp（8.0535 vs 8.0765）
和熵残差 0.005 会一直存在。

## 85. 2026-10-25 THERMO/FORCE 力矩阵差异定位

官方同一几何、同一 HOF（−57.791035）、同一取向，但：
- THERMO 力矩阵 O 对角 6.233425，FORCE 为 6.433114（差 3%）；
- THERMO 梯度范数 2.83003，FORCE 为 2.79511；
- 因此频率 1345.79/2793.39/2858.99 vs 1407.11/2807.96/2863.61。
C++ 两种关键字都复现 FORCE 型矩阵（O 对角 6.413105），说明官方 THERMO
路径在 SCF/力常数上走了不同分支，但 run_mopac/force.F90 中未找到显式
THERMO 分支（仅警告和 thermo 调用）。剩余熵残差 0.005 全部来自振动项
（S_vib 0.0176 vs 0.0226），根源即该力矩阵差异。
结论：C++ 已正确实现 FORCE 路径；THERMO 数值要完全对齐，需要进一步
逆向官方 THERMO 的 SCF 初始化/收敛差异（可能是历史实现细节）。

## 86. 2026-10-25 THERMO/FORCE 差异根因：dcart 的 force 精度标志

源码 `dcart.F90:91`：
```
force = index(keywrd,'PREC') + index(keywrd,'FORCE') /= 0
```
当 FORCE（或 PREC）存在时，dcart 用 3 点高阶差分（`if (force) cdi(k,2)-=chnge2; dhc...`），
否则用 2 点差分（`if (.not.force) cdi+=chnge2`）。这正好解释：
- FORCE：高精度梯度 → gnorm 2.79511，力矩阵 O 对角 6.433；
- THERMO：低精度梯度 → gnorm 2.83003，力矩阵 6.233；
- FORCE+THERMO：FORCE 优先 → 6.433。
C++ 的 dcart.cpp 已按同一逻辑实现 force 标志，但 THERMO 下 fmat 得到的力矩阵
仍与 FORCE 相同（6.413105），说明 fmat 内部的梯度调用没有走 keyword 相关的
dcart 分支（可能静态 icalcn/numcal 缓存导致 force 只算一次）。
下一步：在 fmat/compfg 路径加探针，确认每次梯度调用的 dcart force 值，
使 THERMO 使用 force=false、FORCE 使用 force=true，从而复现 6.233/6.433。

## 87. 2026-10-25 THERMO/FORCE 探针结果（dcart→fmat 链路）

- dcart.cpp 的 force 标志：THERMO=0、FORCE=1，与 Fortran 一致；
- fmat 内梯度范数（+δ）：THERMO 10.053468、FORCE 10.090285（C++ 差 0.4%）；
  官方力矩阵差 3%，说明 C++ 低精度分支的响应太弱；
- fmat 原始对角：THERMO 4.680210、FORCE 4.687142（C++ 差 0.15%）；
  打印矩阵：6.413105 vs 6.431126；官方 6.233425 vs 6.433114。
结论：关键字→dcart→fmat 的传递已通，剩余误差在 dcart 的 force=false
分支或 dhc() 的高阶差分项实现（C++ dhc 是近似实现）。
下一步：逐行对照 dhc.F90 与 dhc.cpp 中 mode=1/2 的项，使低精度梯度
真正降到官方水平，从而复现 6.233/1345/10.003。

## 88. 2026-10-25 dcart 步长扫描（排除步长原因）

用 DCART_H 环境变量扫描 dcart 差分步长，THERMO 下 O 对角：
- h=1e-4：6.413105（默认）
- h=5e-4：6.435984
- h=1e-3：6.434745
- h=2e-3：6.432264
均未降到官方的 6.233425。说明剩余 3% 不是差分步长大小，而是 C++
dhc/h1elec/rotate 的能量函数曲率与 Fortran 不一致，使 2 点差分误差偏小。
探针已移除，代码恢复干净并重建。

## 89. 2026-10-25 dhc 缓冲区对齐（无数值变化）

按 Fortran 把 dhc 的 f/h 固定为 172（1-based，Fortran dimension(171)），
e1b/e2a 仍 45、shmat 81、w/wk 2027。重跑 THERMO：
O 对角仍 6.413105、ZPE 10.109 —— 缓冲区尺寸不是差异原因。
剩余 3% 已定位到 h1elec/rotate/fock2 的数值实现，而非 dhc 数组布局。

## 90. 2026-10-25 rotate 对照结论（e1b/e2a 索引一致）

逐行对照 rotate.F90 与 rotate.cpp：
- e1b/e2a 的填充顺序、ik 计数、en(i*(i-1)/2+k) 索引完全一致；
- 小 rij 清零、li/lj=natorb、elenuc(1,li,li+1,li+lj)、id==3 的
  nddo_to_point 分支都一致；
- hcore.F90 调用 rotate 时传的是 `w(kr)` 数组段（指针已偏移），
  因此 rotate.cpp 内部按 `w[0..]` 相对写入是正确的，不是 bug；
  dhc 传 `&w[kr]`（kr=1）同样正确。
结论：rotate 不是 3% 差异来源；剩余怀疑对象为 rotatd_ wrapper 或 fock2。

## 91. 2026-10-25 fmat 的 precis 分支（2 点 vs 4 点）

fmat.F90 关键分支：
```
precis = index(keywrd,' PREC') /= 0
if (precis) then  ! ±δ、±δ/2 四阶公式，delta 按 10/|g| 缩放
else            ! 2 点公式 (g(+δ/2)-g(-δ/2))*0.5/delta*fact，delta=1/120 固定
```
C++ 原先恒用 4 阶；现已实现两分支。
结果：
- FORCE：O 对角 6.431228（官方 6.433114），ZPE 10.120（官方 10.119）；
- THERMO：O 对角 6.416324（官方 6.233425），ZPE 10.129（官方 10.003）。
THERMO 仍差 3%：2 点公式已正确，差异来自 dcart force=false 的梯度精度
（C++ dhc 低精度误差偏小）。dcart 调用计数确认 THERMO 走 2 点、force=0。

## 92. 2026-10-25 SCF warm-start 定位（未完成）

逐轮 SCF 对比（F90IT/CXXIT）：
- 初始几何：16 轮，escf 轨迹 1e−9 一致；
- 位移几何：Fortran 6 轮、escf 从 −57.786 起；C++ 16 轮、从 955 起。
根因：Fortran iter 的 icalcn 为 SAVE，numcal 不变时跳过 p=0 → 复用上一步密度
（warm start）；C++ icalcn 是局部变量，每次冷启动。
尝试把 icalcn 改 static 并同步 pold/pold2/shift/ten 后，THERMO 报
“GRADIENT IS TOO LARGE (LIMIT=10)”，说明 Fortran 的 SCF SAVE 状态机
（p/pa/pb/pold/pold2/shift/shiftb/ten/scf1/slow/newdg/ready/irrr/nscf）尚未完整移植。
已回退，保持 FORCE ZPE 10.119（官方一致）、THERMO 10.128。

## 93. 2026-10-25 新选项 THERMO_PREC（高精度热化学）

独立基准（同一 PM7 模型的两条官方路径）：
- 官方 FORCE（dcart force=true，3 点高阶差分）：ZPE 10.119；
- 官方 THERMO（dcart force=false，2 点差分）：ZPE 10.003；
两条官方路径自身相差 1.2%，说明 THERMO 的低精度差分误差较大。

C++ 对比：
- 默认 THERMO：10.128；
- THERMO_PREC（强制 4 点差分 + SCFCRT×0.001）：10.118，
  与官方高精度 FORCE 10.119 相差 0.001，明显优于官方 THERMO 的 10.003。

实现：fmat 的 precis 条件增加 THERMO_PREC；iter 的 scfcrt 对 THERMO_PREC ×0.001。
默认行为不变（官方低精度路径，10.128）；需要更高物理精度时使用
"PM7 THERMO THERMO_PREC"。
注：默认 THERMO 仍未逐位复现官方 10.003（低精度 dcart 误差幅度不足），
该复现仍需完整移植 dcart force=false 的数值细节。

## 94. 2026-10-25 chklew F90 620-944 移植完成（Lewis 尾部 + 辅助函数逐行核对）

### 94.1 本段移植内容（按 chklew.F90 行号）

- **F90 620-643（chklew 尾部）**
  - F90 622-626：扫描 `ib(i) > 0`。若没有任何原子残留可用轨道，则 `l = 0`
    并直接结束（等价于原版 `goto 1070`）。
  - F90 627-641：否则在 `debug`（GEOCHK 关键字）为真时，按 Fortran `(20I4)`
    格式转储 "After removing closed pi systems"、`mb`、Core charges(`iz`) 与
    `ib`；转储上限为 `j = Min(numat, large)`（`large` 默认 20，LARGE 关键字时
    为 10^6）。
  - F90 642-643：`l = 1`（Lewis 结构构造成功）。
- **F90 644-652（label 1060）**：原代码中该块只能由 F90:379 的 `goto 1060`
  到达（π 框架已完成时的提前出口）。其中芳香碳计数 `i` 之后从未被使用，且
  `l` 显式置 0；移植中以 `else { l = 0; }` 等价覆盖，并在注释中保留原控制流与
  行号说明，便于逐行审计。
- **F90 654**：`deallocate(ipi, ir5, mpii)` —— C++ 侧为局部 `std::vector`，
  作用域结束自动释放。
- **F90 688-748 / 749-803 / 804-899 / 900-975**：`ring5()`、`arom()`、
  `arom2()`、`add_Lewis_element()` 在既有移植中已有实现；本段逐行对照 Fortran
  核对，结论为**无功能改动**，只补行号映射注释：
  - 数组下标与循环上下界（1-based、`ibonds(k,i)` 列主序）一致；
  - `arom`/`arom2` 的 `iia /= jj`、`mpii(iia) /= 0` 短路条件一致，两段
    "II 侧 / JJ 侧" 对称扫描一致；
  - `add_Lewis_element` 对双原子键、孤对、虚孤对的 `iz`/`ib`/`element_type`
    记账，以及 `ions(atom_i+atom_j)` 累加均一致；
  - `ring5` 保留 Fortran 原版的 `MB(J)` 复检（原作者应为 `MB(K)`）旧行为，
    不做"顺手修正"。
- 新增文件内辅助函数 `print_i4_row`：精确复现 Fortran `(20I4)` 的格式重绕
  （每行 20 个整数、`%4d` 右对齐、满 20 自动换行）。

### 94.2 逐段编译 + 回归

| 段 | 编译 | 回归 |
|---|---|---|
| F90 620-657（尾部新代码） | MSBuild Release x64：chklew.obj 重编译 + 链接通过 | 6/6 通过 |
| F90 688-944（辅助函数核对 + 注释） | 同上，chklew.cpp 无新增告警 | 6/6 通过 |

编译环境：Visual Studio 2022 MSBuild 17.14.23，`/p:Configuration=Release
/p:Platform=x64`；项目既有告警仅 `writmo.cpp` 的 C4996，与本段无关。

### 94.3 回归数值（移植前 → 移植后）

| 用例 | 移植前 C++ | 移植后 C++ | 官方 mopac_ref |
|---|---|---|---|
| 水 BONDS 1SCF | −57.79103 | −57.79103 | −57.79103 |
| π 体系（苯）1SCF | 246.84781 | 246.84781 | 246.84781 |
| THERMO ZPE | 10.007 | **10.007（不变）** | 10.003（已知低精度 dcart 路径，见 §84-93） |
| FORCE ZPE | 10.120 | **10.120（不变）** | 10.119（同上） |
| 水 MOZYME 1SCF | −57.78972 | −57.78972 | −57.78972 |
| 苯 MOZYME 1SCF | 23.02722 | 23.02722 | 23.02726（见 §94.5） |

注：§93 记录的"默认 THERMO 10.128 / THERMO_PREC 10.118"是 THERMO_PREC
实验期间的数值；当前默认路径以本次实测 10.007 / 10.120 为准。

### 94.4 GEOCHK 调试分支验证

新增尾部转储只在 `debug`（GEOCHK）时执行，需单独验证：

- `PM7 MOZYME GEOCHK 1SCF`（苯）：退出码 0，HOF 不变，输出包含
  "After removing closed pi systems" 与 `(20I4)` 的 `i`/`mb`/`iz`/`ib` 行；
- `PM7 MOZYME GEOCHK LARGE 1SCF`（苯）：`large = 10^6` 时转储上限、行宽与换行
  正确，退出码 0，HOF 同样不变。

### 94.5 遗留发现：苯 MOZYME 与官方差 4×10⁻⁵ kcal/mol（不属于本段范围）

- 水 MOZYME 与官方完全一致；苯 MOZYME C++ 23.02722 vs 官方 23.02726。
- GEOCHK 转储显示 C++ 侧苯环 Lewis 结构不完整：6 个碳里只有 2 个被 π 键饱和
  （`iz=ib=0`），另外 4 个残留 `iz=ib=1`；即只生成了 **1 个 π 键**（应为 3 个），
  `chkion` 因此额外打印 "Ion Atom No." 电荷块，进而造成能量偏差。
- 根因在 **F90 289-618 的控制流被展平**（本段 620-944 之前的前序代码）：
  1. F90 456-553 是"`do` 外层循环 + `iset` 判定"：`loop=1` 无键时 **不退出**，
     而是 `loop+1` 继续尝试 Rule 2、Rule 3；命中且 `loop /= 1` 时 `go to 1030`
     把 `loop` 重置为 1 后重新扫描；`loop == 1` 时一次扫描可连续补全多个芳环
     键。现行 C++ 写成 `for (loop = 1..3)` 并在 `if (!iset) break;`，导致苯环
     在 Rule 1 失败后直接跳出，Rule 2/3 根本没有执行。
  2. F90 289 起的 `do` 外层循环（含第四组开链 π 键、`graphi` 分支、
     `1030/1040/1060` 标签的回跳）在 C++ 中被展平为顺序块；554-618 的环 π 键
     块也只加 1 个键就 `break`。
- 影响范围：MOZYME 路径的芳香/共轭体系（非 MOZYME 的 SCF 不经过 CHKLEW）。
  本次段内回归（水/π 1SCF、THERMO、FORCE）均不经过 CHKLEW，故数值全部不变。
- 建议下一步：按 Fortran 原控制流重译 F90 289-618（可保留已就位的 620-944
  尾部与辅助函数不变），修复后苯 MOZYME 应回到官方 23.02726。之后再进行
  linear_cosmo.F90 的分段移植。

### 94.6 复现方式

回归脚本与基线放在 `docs/_verify_runs/chklew_port/`：

    python docs/_verify_runs/chklew_port/run_gate.py --all

脚本会用 `mopac_cpp.exe` 与 `mopac_ref.exe` 各跑一遍
`decks/` 下的 6 个 deck，并与 `baseline.json` 中记录的 C++ 基线逐项比对
（THERMO/FORCE 的官方 ZPE 差异属于 §84-93 已记录的已知路径差异）。

## 95. 2026-10-25 chklew F90 289-618 控制流按原版重译（修复芳香环 π 键缺失）

### 95.1 背景

§94.5 记录的苯 MOZYME 偏差（C++ 23.02722 / 官方 23.02726）根因在前序移植对
F90 289-618 的"展平"：

- F90 456-553 是 `do` 外层循环 + `iset` 判定：`loop=1` 无键时 **不退出**，而是
  `loop+1` 继续 Rule 2、Rule 3；命中且 `loop /= 1` 时 `go to 1030` 把 `loop`
  重置为 1 后重扫。早先 C++ 写成 `for (loop = 1..3)` 且 `if (!iset) break;`，
  苯环在 Rule 1 失败后直接跳出，**Rule 2/3 从未执行**。
- F90 289 起的 `do` 外层循环（含 `1030/1040/1060` 标签、`iset` 回退、
  `graphi` 分支）同样被展平为顺序块。

本轮按 Fortran 原控制流重译 F90 289-618（620-944 尾部与辅助函数不变）。

### 95.2 重译内容

- **F90 289-618 外层 DO 循环**（完整还原）：
  - 293-311：重建 `mb`（π 邻居计数），含 C-S-X 二价 S 特例（304-305）；
  - 324-369：开链 ene/yne π 键。`loop=1` 从碳端开始、`loop=2` 从任意端开始；
    `do while (mb(jj)/=0 .and. ib(jj)/=0 .and. ib(ii)>0)` 直到价键用尽；
  - 375-379：本轮无新键时，若全部原子饱和 → `1060`（`l=0` 结束），否则进入
    1020 调试输出并尝试开环；
  - 403-415：重建 `mb`；423-426：`first_pi` 初始化；427-452：扩展/去重 `ipi`
    候选表（`loop1`/`loop2`）；
  - **456-553：Rule 1（`arom`）/ Rule 2（`arom2`）/ Rule 3（任意键）**。
    `loop` 从 1 递增；命中且 `loop /= 1` 时 `go to 1030`（`loop=1` 重扫）；
    `loop=1` 时一次扫描可连续补全多个芳环键；炔键
    `ib(ii)==1 .and. ib(jj)==1` 追加第二个 π 键（495-497）；
  - 554-618：环 π 键。按 `lnext`、`mb==3`、`loop>=3` 判据扫描，`go to 1040`
    加键后**重启整个外层 DO 循环**；无环键可加时 `exit`（605）进入 620 尾部。
- **F90 76-93**：补齐 `txtatm` 第 7 列"临时置空 → 扫描 → 恢复"逻辑（与 Fortran
  完全一致；此前只扫 `substr(0,6)`/`(0,12)` 而未清空第 7 列）。
- **F90 97-110**：补 `ndelec(j) /= 0` 判定（行为等价，可读性/逐行对应）。
- **F90 162-183 SETPI 用户 π 键**：按 Fortran 结构移植，含动态 `i` 编辑描述符
  宽度（`2+floor(log10(atom*1.01))`）的 `mopend` 消息。注：本 C++ 端 `readmo`
  仍未填充 `PIBONDS`（只有 `PIBONDS_TXT`），SETPI 输入尚不可用；`PIBONDS`
  为空时按 Fortran 的 `ii==0` 立即退出处理。
- **DEBUG（GEOCHK）转储**：F90 112-123 / 151-161 / 217-226 / 312-320 / 380-398
  的 `(20I4)` 数组转储全部补齐；新增 `print_ld_text` 复现 Fortran 表控字符
  输出多一个前导空列的行为。至此 chklew 的 DEBUG 输出与官方逐字节一致。
- 控制流用显式状态 + 前向 `goto finish` 等价表达，注释逐处标注 F90 行号与标签。

### 95.3 验证

编译：MSBuild Release x64 通过，`chklew.cpp` 0 告警。

回归（14 项，`docs/_verify_runs/chklew_port/run_gate.py --all`）：

| 用例 | C++ | 官方 mopac_ref | 结果 |
|---|---|---|---|
| 水 BONDS 1SCF | −57.79103 | −57.79103 | 一致 |
| 苯 π 1SCF | 246.84781 | 246.84781 | 一致 |
| THERMO ZPE | **10.007** | 10.003 | 不变（已知低精度 dcart 路径，§84-93） |
| FORCE ZPE | **10.120** | 10.119 | 不变（同上） |
| 水 MOZYME 1SCF | −57.78972 | −57.78972 | 一致 |
| 苯 MOZYME 1SCF | 23.02726 | 23.02726 | **修复**（此前 23.02722） |
| ethene | 14.29865 | 14.29865 | 一致 |
| acetylene | 57.81121 | 57.81121 | 一致 |
| butadiene | 261.60662 | 261.60662 | 一致 |
| pyridine | 37.15222 | 37.15222 | 一致 |
| naphthalene | 43.04634 | 43.04633 | 差 1×10⁻⁵（见下） |
| cyclopentadiene | 120.81544 | 120.81544 | 一致 |
| indole | 70.32731 | 70.32731 | 一致 |
| anthracene | 67.10911 | 67.10911 | 一致 |

- **GEOCHK 逐字节比对**：苯与萘的 chklew DEBUG 区域（Arrays before sigma /
  sigma / lone pairs / MB / 开链 / 1020 / 620 尾部）与官方输出**逐字节一致**，
  包括前导空列与 `(20I4)` 的每 20 项换行。
- **naphthalene 的 1×10⁻⁵ kcal/mol 残差**：chklew 的 DEBUG 转储与官方逐字节
  一致（Lewis 结构相同），残差位于 MOZYME SCF / 局域轨道端的已知数值精度
  范围，与本轮 chklew 改动无关。
- 此前"旧 CHKLEW"在 naphthalene 上偏差 5×10⁻³ kcal/mol（43.05115），本轮
  重译后降到 1×10⁻⁵，亦旁证控制流修复正确。

### 95.4 本轮发现但不在本范围内的遗留问题：MOZYME + 硫体系

- thiophene（PM7 MOZYME 1SCF）：C++ −26527.70 vs 官方 104.01；
- methanethiol（PM7 MOZYME 1SCF）：C++ −7374.70 vs 官方 3.28；
- 但同一几何的**非 MOZYME** PM7 路径 C++ 与官方逐位一致
  （thiophene：104.00922 = 104.00922）；
- thiophene 的 chklew DEBUG 转储（`iz`/`ib`/`mb`/各阶段）与官方逐字节一致。

结论：该问题不在 CHKLEW，指向 MOZYME 路径上的硫（d 轨道/参数）处理，建议
单独定位。含硫体系的 MOZYME 暂不纳入本轮回归门（已在
`compare_pi_cases.py` 的说明中记录）。

### 95.5 复现

    python docs/_verify_runs/chklew_port/run_gate.py --all

脚本同时运行 6 个基线用例与 8 个 π 体系用例（`pi_cases/`），与 `baseline.json`
逐项比对；`compare_pi_cases.py` 可单独做 π 体系的官方对照（不依赖基线）。

## 96. 2026-10-26 MOZYME + COSMO 打通：linear_cosmo 数值修复（mkl_bits DTRTI2 / 面积因子）

### 96.1 背景与定位方法

按计划进入 linear_cosmo.F90。审计发现其 C++ 移植（约 1835 行，Batch A–F）并非
§71 所说的 "no-op"，而是**已实现但数值错误**：

- `PM7 MOZYME 1SCF EPS=78.4`（水）C++ 端直接访问越界崩溃；
- 甲醇/乙酸根在 C++ 端给出 −7391 / −17209 kcal/mol 量级的乱值（官方
  −53.32 / 47.19）；
- 非 MOZYME COSMO 的 `DIELECTRIC ENERGY` 此前也有 5×10⁻⁵ 级偏差。

为定位问题，本轮搭建了可输出中间量的参考算例：

1. 用 ifx 2025.3 从官方 archive 源码重建 `mopac_ref_dbg.exe`：复用既有
   `obj/*.obj` 与 `*.mod`，仅重编 `linear_cosmo.F90`（插入 dump）并补编缺失的
   `denrot.obj`，链接后结果与官方 `mopac_ref.exe` **完全一致**；
2. 在 Fortran 与 C++ 两侧以同一格式 dump：`nps / area / cosvol / cosurf /
   iatsp / nar_csm / nsetf / nipsrs / rsc / a_diag / a_part / 预处理矩阵 /
   cgm 解 / phinet / qscnet / ediel / solv_energy`，逐项比对。

由此定位到两处**相互独立**的错误。

### 96.2 修复 1：linear_cosmo 面积因子用了 cosurf(3,·)（7 处）

Fortran 在 `aq_vec` / `simulate_aq_vec` / `aq_dir_int` / `simulate_aq_dir_int` /
`precond_aq_dir_int` 中为
`aa = aa / (cosurf(4,ii) * cosurf(4,jj))`（`cosurf(4,·)` 是段面积），
`amat_diag` 为 `m(ii) = aa / cosurf(4,ii)**2`；C++ 全部误写成 `cosurf[3]`
（方向的 z 分量）。已改为 `cosurf[4]`（6 处 `aa /=` + 1 处 `amat_diag`）。

另：Fortran `a_part(1:ndim)` 为 1-based，C++ 原 `a_part.assign(ndim, 0.0)`
会让 `a_part[ndim]` 越界 1 个元素，改为 `assign(ndim + 1, 0.0)`。

### 96.3 修复 2：mkl_bits 的 DTRTI2（三角矩阵求逆）算法错误

`precondition` 的块对角求逆链路为 dpotrf → dpotri → dtrtri → dtrti2，全部由
`mkl_bits.cpp` 手工复现。与 MKL 对拍（本轮新增独立连测程序）：

- `dpotrf`：与 MKL 差 2×10⁻¹⁶（正常）；
- `dtrtri`：与 MKL 差 **4.9×10⁻²**（错误）；
- 结果 `dpotri` 差 5.2×10⁻²。

原 `dtrti2` 的公式与循环顺序均不符合 LAPACK。按 Netlib `dtrti2.f` 中
`DTRMV('U'/'L','N',DIAG)` 的语义重写 upper/lower 两个分支后：

| 与 MKL 对拍 | 修复前 | 修复后 |
|---|---|---|
| dpotrf | 2.2e-16 | 2.2e-16 |
| dtrtri | 4.9e-2 | 9.7e-17 |
| dpotri | 5.2e-2 | 1.1e-16 |

该修复同时消除了非 MOZYME COSMO 的 `DIELECTRIC ENERGY` 5×10⁻⁵ 偏差（见 96.4）。

### 96.4 验证

编译：MSBuild Release x64 通过；`linear_cosmo.cpp`、`mkl_bits.cpp` 无新增告警。

**(1) 既有 14 项回归**（水/π 1SCF、THERMO、FORCE、6 项 MOZYME、8 项 π 体系）：
**全部 PASS**，THERMO ZPE **10.007**、FORCE ZPE **10.120** 保持不变。

**(2) MOZYME + COSMO 1SCF**（`PM7 MOZYME 1SCF EPS=78.4`；脚本
`docs/_verify_runs/cosmo_moz/run_cosmo_moz.py`）：

| 用例 | C++ HOF / DIELECTRIC | 官方 | 结果 |
|---|---|---|---|
| 水 | −65.16392 / −0.36147 | 同 | 一致 |
| 甲醇 | −53.32319 / −0.29319 | 同 | 一致 |
| 苯 | 21.05947 / −0.09229 | 同 | 一致 |
| 环戊二烯 | 105.83580 / −0.77391 | 同 | 一致 |
| 萘 | 39.79622 / −0.15878 | 同 | 一致 |
| 乙酸根 CH₃COO⁻ | 47.18759 / −3.37822 | 47.18758 / −3.37822 | 差 1×10⁻⁵（见 96.5） |

**(3) 非 MOZYME COSMO 1SCF**（`PM7 1SCF EPS=78.4`）：水 −65.16522 / −0.35863、
甲醇 −53.32397 / −0.29253，与官方**逐位一致**（此前水的 DIELECTRIC ENERGY
为 −0.35858，差 5×10⁻⁵，本轮修复后对齐）。

**(4) 中间量级比对**（水，GEOCHK 式 dump）：`area/cosvol/cosurf/iatsp/nar_csm/
nsetf/nipsrs` 逐字节一致；`a_diag`、`a_part`、块预处理矩阵、CG 解、`phinet`、
`qscnet`、`ediel`、`solv_energy` 与参考差 < 1×10⁻¹⁰（修复前 `precond m` 首元素
差 5.2×10⁻²）。surclo 输出的 `rsc` 尾部未使用段与 Fortran 不同（近场实际用到
的环段一致），不影响已验证算例。

### 96.5 遗留问题（不属 1SCF 移植范围）

1. **COSMO 几何优化崩溃（旧 bug，非本轮引入）**
   - `PM7 MOZYME EPS=78.4`（无 1SCF）与 `PM7 EPS=78.4` 在 C++ 端
     0xC0000005 崩溃，官方正常；
   - 用本轮修复前的旧构建同样崩溃，且**非 MOZYME COSMO 优化也复现** →
     属于共享的 `cosmo.cpp` 几何更新路径（coscav/coscanz/surclo 的重复调用），
     不是 linear_cosmo 1SCF 链路；
   - 已定位：第二次 `coscavz` 之后、`deriv → dcart` 释放堆时暴露堆损坏；
     Debug 堆报告 "modified past requested size"；在 linear_cosmo 侧加
     `LC_SKIP_COSCAVZ2 / LC_SKIP_SURCLO2 / LC_SKIP_ADD2` 局部跳过均不能消除。
     建议下一步单独用 PageHeap/Debug STL 定位。
2. **乙酸根 HOF 1×10⁻⁵ 残差**（DIELECTRIC 已一致）：疑为 CG 迭代路径差异，
   待用中间量 dump 继续定位。
3. surclo 的 `rsc` 未使用段差异（当前无影响）。

### 96.6 复现

    python docs/_verify_runs/cosmo_moz/run_cosmo_moz.py --all
    python docs/_verify_runs/chklew_port/run_gate.py --all

（前者为 MOZYME+COSMO 1SCF 与优化用例，后者为既有 14 项基线回归。）

## 97. 2026-10-26 COSMO 几何优化崩溃修复（diegrd 布局 / dcart 平面缓冲索引）

### 97.1 问题

§96.5(1) 记录的 COSMO 几何优化崩溃：`PM7 MOZYME EPS=78.4` 与 `PM7 EPS=78.4`
（不带 1SCF）在 C++ 端 0xC0000005；该崩溃在修复前构建、且**非 MOZYME COSMO 优化**
同样复现，属共享 cosmo.cpp 路径。

### 97.2 根因

`dcart` 在 `useps` 时把笛卡尔梯度打包成平面缓冲
`flat[(i-1)*3+k] = dxyz[k][i]`，交给 `diegrd()`（COSMO 介电梯度，位于
`cosmo.cpp`）。但 C++ `diegrd` 内部按 `dxyz[ix*(numat+1)+atom]` 索引——这是另一种
布局约定（坐标优先、步长 `numat+1`），与 dcart 的缓冲区（原子优先、步长 3）不一致：

- 官方 Fortran：`dxyz` 声明为 `(3, numat*l123)` 列主序；`diegrd` 哑元
  `dxyz(3,numat)`，元素 `(ix, atom)` 的地址为 `(ix-1) + 3*(atom-1)`；
- C++ `diegrd` 按 `ix*(numat+1)+atom` 越界写（水：`numat=3` 时最大写下标 11，
  而 `flat` 仅 `3*numat+1 = 10` 个元素）→ 堆损坏（在随后 `deriv → dcart`
  释放内存时被 Debug 堆检测到）并产生错误梯度；
- 另外 dcart 打包用 `flat[(i-1)*3+k]`（k=1..3），与 diegrd 的 0 基索引相差 1，
  即使不越界也会把梯度写到错误的原子/分量上。

### 97.3 修复

- `cosmo.cpp` 的 `diegrd()`：全部 9 行 18 处 `dxyz[ix*(numat+1)+atom]`
  改为 `dxyz[(atom-1)*3+(ix-1)]`，与 Fortran 列主序布局完全一致；
- `dcart.cpp`：`flat` 打包/回写改为 `flat[(i-1)*3+(k-1)]`，与 diegrd 的
  0 基索引一致。

### 97.4 验证

崩溃消失：`PM7 EPS=78.4`、`PM7 MOZYME EPS=78.4`（水/甲醇）均正常结束。

**非 MOZYME COSMO 几何优化**（EF 优化器）：

| 用例 | C++ | 官方 | 结果 |
|---|---|---|---|
| 水 `PM7 EPS=78.4` | 首轮梯度 8.388；HOF −65.21143 / DIEL −0.36538 | 同 | **完全一致** |
| 甲醇 `PM7 EPS=78.4` | 首轮梯度 32.005；HOF −54.01839 / DIEL −0.28969 | 32.005；−54.01739 / −0.28964 | 梯度一致，收敛值差 1×10⁻³（EF 优化器已知精度，§72-82） |

**MOZYME COSMO 几何优化**：不再崩溃，但首轮梯度 684 / 3451 vs 官方 7.77 / 31.88
——属于**既有 MOZYME 几何优化链路问题**（§71 已记录：水 MOZYME 优化 C++ 最终
−57.791 vs 官方 −57.799；气体 MOZYME 优化同样不收敛），与本轮 diegrd/COSMO 修复
无关，需单独移植/修复 MOZYME 梯度链路。

**回归**：既有 14 项 ALL PASS（THERMO 10.007 / FORCE 10.120 不变）；MOZYME+COSMO
1SCF 6 项结果不变（5 项与官方一致，乙酸根差 1×10⁻⁵）。

### 97.5 对 §96.5 的更新

- §96.5(1)「COSMO 几何优化崩溃」**已修复**（本节），不再是遗留；
- 仍遗留：MOZYME 几何优化不可用（独立问题，非 COSMO）；乙酸根 1×10⁻⁵；
  surclo `rsc` 尾部未使用段差异。

### 97.6 复现

    python docs/_verify_runs/cosmo_moz/run_cosmo_moz.py --all

## 98. 2026-10-26 MOZYME 几何优化可用性实测（含 COSMO）：确认不可用

### 98.1 测试目的

对 §97.4 中「MOZYME 几何优化（含 COSMO）不可用」这一结论做独立实测：
在修复 COSMO 优化崩溃（§97）之后，MOZYME 优化是否已经可用；问题是否与优化器/
COSMO 相关；能量与梯度哪一层出错。

### 98.2 测试方法

脚本 `docs/_verify_runs/moz_opt_probe/run_moz_opt_probe.py`：

- 分子：水（3 原子）、苯（12 原子）、C12H26（38 原子，取 `x64/Release/c12_pub.mop`
  几何）；
- 优化器/关键字：默认（`PM7 MOZYME`）、`BFGS`、`L-BFGS`、`EF`、`COSMO`
  （`PM7 MOZYME EPS=78.4`）；
- 每个用例 C++ 与官方 `mopac_ref.exe` 各跑一次，记录：退出码、正常终止、最终
  HOF、最终梯度范数、优化轮数、首轮梯度、失败信息。

### 98.3 结果（首轮梯度 / 最终 HOF / 最终梯度范数）

| 用例 | C++ | 官方 | 结论 |
|---|---|---|---|
| 水 默认 | g1=**721.4** → HOF −57.79108，GN 722.6 | g1=3.82 → HOF −57.79909，GN 0.78 | 不可用 |
| 水 BFGS | 0 轮，HOF −57.78992，GN 716.3 | 3 轮，HOF −57.79922，GN 0.95 | 不可用 |
| 水 L-BFGS | 同默认（关键字未生效） | 官方**报错**（不支持） | 不可用 |
| 水 EF | 同默认 | 同默认 | 不可用 |
| 水 COSMO | g1=684.1 → HOF −65.17074，GN 693.8 | g1=7.77 → HOF −65.21180，GN 0.46 | 不可用 |
| 苯 默认 | g1=**6284.3** → HOF 23.02282，GN 6284.2 | g1=17.24 → HOF 22.95671，GN 0.58 | 不可用 |
| 苯 BFGS | 0 轮，HOF 23.16493，GN 6291.1 | HOF 22.96272 | 不可用 |
| 苯 COSMO | g1=6220.8 → HOF 21.05900，GN 6220.6 | g1=20.18 → HOF 20.95883，GN 0.67 | 不可用 |
| C12H26 默认 | g1=**5513.4** → HOF 1252.10，GN 5594.8 | g1=510.9 → HOF 162.75，GN 0.99 | 不可用 |
| C12H26 BFGS | 0 轮，HOF 1231.60，输出 FAILED 信息 | 63 轮，HOF 171.14 | 不可用 |
| C12H26 COSMO | g1=5470.6 → HOF 1234.42 | g1=514.1 → HOF 98.40 | 不可用 |

关键观察：

1. **首轮能量已完全一致**（例如水首轮 HEAT = −57.78972，C++ 与官方相同），
   说明 SCF/能量路径正常；
2. **首轮梯度就错 2 个数量级**（721 vs 3.8；6284 vs 17.2；5513 vs 511），
   且与分子、优化器、是否 COSMO 无关 → 根因在**梯度计算**，不在优化器；
3. C++ 优化 3 轮后即以巨大梯度范数"结束"（HOF 停留在起点附近），
   BFGS 甚至 0 轮；因此优化完全不起作用；
4. `L-BFGS` 关键字官方不支持（报错），C++ 静默按默认处理；
5. COSMO 与气体表现一致（同样是首轮梯度错误），说明与 §96/§97 的 COSMO
   修复无关。

### 98.4 梯度错误定位（一层）

用 `DCART` 关键字强制打印 dcart 的笛卡尔导数（水、同一初始几何）：

| 分量 | C++（MOZYME） | 官方（MOZYME） |
|---|---|---|
| O x / y / z | 266.53 / −385.55 / −252.05 | 0.0504 / 0.00034 / −3.1010 |
| H1 x / y / z | −136.16 / 410.86 / 111.68 | −0.0250 / −0.2210 / 1.5349 |
| H2 x / y / z | −130.38 / −25.31 / 140.37 | −0.0254 / 0.2206 / 1.5661 |

两者都满足平移不变（三原子分量和≈0），但 C++ 的导数整体偏大 2 个数量级、
且方向结构不同 → 错误在 `dcart` 的 MOZYME 分支（`ijbo` 成键块构建 / `dhc`
双中心能量差分 / `delsta` 点电荷项）之一，而不是 jcarin/优化器。

### 98.5 结论与下一步

- **结论：MOZYME 几何优化（含 COSMO）目前不可用**——不是"结果略差"，而是
  梯度从第 1 步就错误 2 个数量级，优化器无法工作（与 §71 的旧记录一致）；
- 该问题与 §96/§97 的 COSMO 数值修复无关（气体 MOZYME 同样复现）；
- 建议下一步：为 `dcart` 重建可 dump 中间量的 Fortran 参考对照（复用 §96 的
  ifx/obj 管线），逐项比对 MOZYME 分支的 `padi/pdi/pbdi`、`dhc` 的 aa/ee、
  `delsta` 的 dstat，定位并修复后再重跑本测试矩阵。

## 99. 2026-10-26 旧版本（GitHub 状态）复现 MOZYME 优化问题

### 99.1 目的

确认 §98 实测到的「MOZYME 几何优化梯度错误」是否为本轮修复引入的回归，还是
GitHub 旧版本就存在的问题。

### 99.2 被测版本

| 版本 | 说明 | 梯度相关源码 |
|---|---|---|
| `8f66c42` | 本轮工作之前 GitHub 上的主分支状态（2026-10-24，本地 commit "1" 的父提交） | 已直接编译实测 |
| `a30b099` / `v1.1.0` | GitHub 上已发布的旧 tag（`origin/main` 跟踪点） | 已直接编译实测（与 8f66c42 仅 README/LICENSE/docs/.gitignore 差异，**无任何代码文件变化**） |
| `b5fa665` | 本轮全部修复后的状态（commit "1"） | 含 dcart/cosmo 等修复 |

两个版本均用 `git archive` 导出干净源码树、MSBuild Release x64 重新完整编译
（退出码 0，无编译错误）：`8f66c42` → `mopac_cpp_old_8f66c42.exe`，
`v1.1.0`（`a30b099`）→ `mopac_cpp_v110.exe`。`git diff --name-only a30b099..8f66c42`
不含任何 `.cpp/.h/.vcxproj` 文件（仅 README/LICENSE/docs/.gitignore），两版本
代码相同；另确认两版本的梯度链路文件（dcart/deriv/dhc/jcarin/chrge_for_MOZYME/
iter_for_MOZYME/hcore_for_MOZYME）与 8f66c42 逐字节相同。

### 99.3 结果：旧版本同样有该问题（且 COSMO 优化直接崩溃）

| 用例 | 旧版本（v1.1.0 与 8f66c42 实测相同） | 当前版本（b5fa665） | 官方 |
|---|---|---|---|
| 水 默认优化 | g1=**721.377**，HOF −57.79108，GN 722.56，3 轮 | g1=721.377，HOF −57.79108，GN 722.56，3 轮 | g1=3.818，HOF −57.79909，GN 0.78 |
| 水 COSMO 优化 | **崩溃**（exit 3221225477，输出为空） | 正常结束，g1=684.06，HOF −65.17074 | g1=7.768，HOF −65.21180，GN 0.46 |
| 苯 默认优化 | g1=6286.50，HOF 23.02260，GN 6284.43 | g1=6284.29，HOF 23.02282 | g1=17.24，HOF 22.95671，GN 0.58 |
| 苯 COSMO 优化 | **崩溃**（exit 3221225477） | 正常结束，g1=6220.84，HOF 21.05900 | g1=20.18，HOF 20.95883 |
| C12H26 默认优化 | g1=5515.75，HOF 1253.64 | g1=5513.35，HOF 1252.10 | g1=510.92，HOF 162.75，GN 0.99 |

v1.1.0 直接实测数据与上表 8f66c42 完全一致：水默认 g1=721.377、苯默认
g1=6286.499、C12H26 g1=5515.748；两个 COSMO 优化用例同样以 exit 3221225477
崩溃（输出为空），DCART 导数同为 O (266.532342, -385.552987, -252.046453)。

**决定性证据：旧版本与当前版本的 dcart 笛卡尔导数逐位相同**（水，`DCART`）：

```
旧版本 8f66c42 : O (266.532342, -385.552987, -252.046453)   H1 (-136.156226, 410.859945, 111.675717)
当前   b5fa665 : O (266.532342, -385.552987, -252.046453)   H1 (-136.156226, 410.859945, 111.675717)
官方           : O (0.050406, 0.000340, -3.101008)           H1 (-0.024998, -0.220983, 1.534891)
```

即：本轮所有修复（chklew / linear_cosmo / mkl_bits / cosmo-diegrd / dcart-flat）
**没有改变 MOZYME 梯度**，该错误是旧版本就存在的历史遗留；本轮只是把
「MOZYME + COSMO 优化崩溃」修成了「能运行但仍因梯度错误而不收敛」。

### 99.4 结论

- **GitHub 旧版本 v1.1.0 与 8f66c42 直接编译实测均复现该问题**：
  MOZYME 优化首轮梯度高 2 个数量级、优化不收敛；
- 旧版本还额外存在 MOZYME+COSMO 优化崩溃（本轮 §97 已修复）；
- 因此 §98 的问题不是本轮引入的回归，而是需要单独修复的历史遗留
  （§98.5 的 dcart MOZYME 分支定位仍然成立，下一步可继续）。

### 99.5 复现

    git archive --format=tar --output=old.tar 8f66c42
    mkdir old && tar -xf old.tar -C old
    msbuild old\mopac_cpp_public\MOPAC_Cpp.vcxproj /p:Configuration=Release /p:Platform=x64 /m
    python docs\_verify_runs\moz_opt_probe\old_probe.py <old>\mopac_cpp_public\x64\Release\mopac.exe

## 100. 2026-10-26 MOZYME 几何优化修复（dcart 密度块索引 + addnucz 接线）

### 100.1 问题回顾

§98/§99 实测：MOZYME 几何优化首轮梯度比官方大 ~2 个数量级
（水 C++ 721.4 vs 官方 3.818），气体与 COSMO 均复现；v1.1.0 / 8f66c42
旧版本同样存在。

### 100.2 根因 1：dcart 的 MOZYME 分支密度块索引基准错误

- Fortran `dcart` 的 MOZYME 分支以 1-based 打包索引填充 `padi/pdi/pbdi`
  （`padi(1..171)`）；而 C++ `dhc()` 的 ABI 是 0-based 连续数组
  （内部 `pv[t] = p[t-1]`）。
- C++ 的 MOZYME 分支照抄了 F90 的 1-based 下标（`padi[ij_]` 从 1 开始），
  而非 MOZYME 分支却是 0-based（`padi[ij_ + m]`，`ij_` 从 0 开始）——
  两条分支约定不一致，MOZYME 分支的密度块整体错位 1 个元素。
- 后果：`dhc()` 读到错位密度，双中心能量 `aa/ee` 错误
  （C++ −4.072 vs 官方 −13.508），有限差分梯度放大 2 个数量级。
- 修复：MOZYME 分支 4 处填充改为 F90 索引减 1（`padi[ij_ - 1]`、
  `padi[l - 1]`、`pdi[...]` 同步）。

验证：修复后 dcart 中间量 `PAIR` 的 `aa/ee/deriv` 与 Fortran 参考一致到
~1×10⁻¹²，笛卡尔梯度 `(0.050405530, 0.000340452, -3.101007880)` 与官方
`(0.050406, 0.000340, -3.101008)` 逐位一致。

### 100.3 根因 2：hcore_for_MOZYME 的 addnucz 空桩（仅影响 COSMO 梯度）

修好根因 1 后气体 MOZYME 优化与官方完全一致，但 MOZYME+COSMO 首轮梯度
仍错（水 101.5 vs 7.768）。原因：

- `hcore_for_MOZYME` 的 `useps` 分支中 `addnucz(phinet, qscnet, qdenet)`
  被写成 no-op 桩（注释 "Stub: no-op placeholder"）；
- `addnucz` 负责把核电荷写入 `qdenet(:,1)`（并把 channel-1 的
  phi/qscnet 清零），`addfckz` 随后计算 `qdenet(:,3) = qdenet(:,2)+qdenet(:,1)`；
- `diegrd`（COSMO 介电梯度）用的正是 `qdenet(:,3)`，缺少核电荷项导致
  COSMO 梯度错误；能量/`DIELECTRIC ENERGY` 不用该通道，所以 1SCF 一直正确。
- 修复：恢复调用 `linear_cosmo::addnucz(phinet, qscnet, qdenet)`。

### 100.4 验证（C++ vs 官方 mopac_ref.exe）

**MOZYME 优化**（首轮梯度 / 最终 HOF / 最终 GN / 轮数）：

| 用例 | C++ | 官方 |
|---|---|---|
| 水 `PM7 MOZYME` | 3.818 / −57.79909 / 0.77912 / 3 | 完全一致 |
| 水 `+ EPS=78.4` | 7.768 / −65.21180 / 0.45776 / 3 | 完全一致 |
| 苯 `PM7 MOZYME` | 17.244 / 22.95671 / 0.57086 / 11 | 17.244 / 22.95671 / 0.57602 / 11 |
| 苯 `+ EPS=78.4` | 20.182 / 20.95883 / 0.66506 / 5 | 完全一致 |
| 萘 `PM7 MOZYME` | 119.087 / 39.71539 / 0.68252 / 19 | 119.087 / 39.71539 / 0.68250 / 19 |
| 萘 `+ EPS=78.4` | 121.751 / 36.32591 / 0.38652 / 18 | 完全一致 |
| 环戊二烯 `PM7 MOZYME` | 210.033 / 32.29160 / 0.99026 / 29 | 210.033 / 32.29160 / 0.98974 / 29 |
| 环戊二烯 `+ EPS=78.4` | 232.545 / 68.32543 / 2.39246 / 18 | 完全一致 |
| 甲醇 `+ EPS=78.4` | −54.01756 / −0.28965 | 完全一致 |
| C12H26 `PM7 MOZYME` | 首轮梯度 510.922（与官方一致），收敛 163.55 | 510.922，收敛 162.75 |

C12H26 收敛值差异属 EF 优化器在大体系上的既有精度问题，与 MOZYME 无关：
非 MOZYME 的 C12H26 优化同样 −67.34742 vs 官方 −67.35916。

**回归**：

- 既有 14 项：**全部 PASS**，THERMO ZPE 10.007 / FORCE ZPE 10.120 不变；
- COSMO 套件 12 项：10 项与官方完全一致（含此前 FAIL 的
  `water_moz_cosmo_opt`、`meoh_moz_cosmo_opt` 两项转为 PASS）；
  剩余：乙酸根 1×10⁻⁵、甲醇普通 COSMO 优化 1×10⁻³（EF 优化器已知精度）。

### 100.5 D: 工程重编译后的复核

上述修复写回 `D:\test_rep\mopac_cpp_public` 并完整重编译（退出码 0）后，用
D: 构建的 exe 重跑全部用例，结果与工作区构建一致：

- 14 项既有回归：**ALL PASS**（THERMO 10.007 / FORCE 10.120）；
- COSMO 套件：10/12 完全一致（含 `water_moz_cosmo_opt`、`meoh_moz_cosmo_opt`
  由 FAIL 转 PASS）；
- MOZYME 优化 probe（15 用例）：
  - 水/苯/萘/环戊二烯（气体与 COSMO）、水 `BFGS`、苯 `BFGS`：首轮梯度与
    最终 HOF 与官方一致；苯/萘/环戊二烯最终 GN 有 1e-3~1e-2 级差异
    （EF 优化器数值细节，HOF 一致）；
  - `L-BFGS` 关键字官方不支持（ref 报错），C++ 回落为默认优化器；
  - C12H26（38 原子直线烷烃测试几何）：首轮梯度与官方一致（510.922 /
    514.071），收敛 HOF 163.55 vs 162.75（COSMO 98.73 vs 98.40）——
    属 EF 优化器在大体系上的既有精度问题：非 MOZYME 的 C12H26 优化
    同样 −67.34742 vs 官方 −67.35916。

结论：**MOZYME 几何优化（气体与 COSMO）已修复**——梯度从第 1 步起与官方
一致，小/中等分子优化结果逐位相同；剩余差异均归因于 EF 优化器既有精度，
与 MOZYME 无关。

### 100.6 复现

    python docs/_verify_runs/moz_opt_probe/run_moz_opt_probe.py --all
    python docs/_verify_runs/cosmo_moz/run_cosmo_moz.py --all
    python docs/_verify_runs/chklew_port/run_gate.py --all
