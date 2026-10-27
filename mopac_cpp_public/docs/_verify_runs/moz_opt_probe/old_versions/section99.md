## 99. 2026-10-26 旧版本（GitHub 状态）复现 MOZYME 优化问题

### 99.1 目的

确认 §98 实测到的「MOZYME 几何优化梯度错误」是否为本轮修复引入的回归，还是
GitHub 旧版本就存在的问题。

### 99.2 被测版本

| 版本 | 说明 | 梯度相关源码 |
|---|---|---|
| `8f66c42` | 本轮工作之前 GitHub 上的主分支状态（2026-10-24，本地 commit "1" 的父提交） | — |
| `a30b099` / `v1.1.0` | GitHub 上已发布的旧 tag（`origin/main` 跟踪点） | 与 `8f66c42` **逐字节相同** |
| `b5fa665` | 本轮全部修复后的状态（commit "1"） | 含 dcart/cosmo 等修复 |

用 `git archive 8f66c42` 导出干净源码树、MSBuild Release x64 重新完整编译
（退出码 0，无编译错误），得到 `mopac_cpp_old_8f66c42.exe`；
`git diff a30b099..8f66c42 -- dcart.cpp deriv.cpp dhc.cpp jcarin.cpp
chrge_for_MOZYME.cpp iter_for_MOZYME.cpp hcore_for_MOZYME.cpp` 为空，
即 `v1.1.0` 的梯度链路与本测试版本完全相同。

### 99.3 结果：旧版本同样有该问题（且 COSMO 优化直接崩溃）

| 用例 | 旧版本 8f66c42 | 当前版本（b5fa665） | 官方 |
|---|---|---|---|
| 水 默认优化 | g1=**721.377**，HOF −57.79108，GN 722.56，3 轮 | g1=721.377，HOF −57.79108，GN 722.56，3 轮 | g1=3.818，HOF −57.79909，GN 0.78 |
| 水 COSMO 优化 | **崩溃**（exit 3221225477，输出为空） | 正常结束，g1=684.06，HOF −65.17074 | g1=7.768，HOF −65.21180，GN 0.46 |
| 苯 默认优化 | g1=6286.50，HOF 23.02260，GN 6284.43 | g1=6284.29，HOF 23.02282 | g1=17.24，HOF 22.95671，GN 0.58 |
| 苯 COSMO 优化 | **崩溃**（exit 3221225477） | 正常结束，g1=6220.84，HOF 21.05900 | g1=20.18，HOF 20.95883 |
| C12H26 默认优化 | g1=5515.75，HOF 1253.64 | g1=5513.35，HOF 1252.10 | g1=510.92，HOF 162.75，GN 0.99 |

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

- **GitHub 旧版本（8f66c42；v1.1.0 梯度源码与之相同）同样存在该问题**：
  MOZYME 优化首轮梯度高 2 个数量级、优化不收敛；
- 旧版本还额外存在 MOZYME+COSMO 优化崩溃（本轮 §97 已修复）；
- 因此 §98 的问题不是本轮引入的回归，而是需要单独修复的历史遗留
  （§98.5 的 dcart MOZYME 分支定位仍然成立，下一步可继续）。

### 99.5 复现

    git archive --format=tar --output=old.tar 8f66c42
    mkdir old && tar -xf old.tar -C old
    msbuild old\mopac_cpp_public\MOPAC_Cpp.vcxproj /p:Configuration=Release /p:Platform=x64 /m
    python docs\_verify_runs\moz_opt_probe\old_probe.py <old>\mopac_cpp_public\x64\Release\mopac.exe
