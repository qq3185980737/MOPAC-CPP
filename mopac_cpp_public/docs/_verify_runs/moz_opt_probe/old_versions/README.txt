old_versions — 旧版本（v1.1.0 / 8f66c42）MOZYME 优化复现测试（2026-10-26）
=====================================================================

直接编译实测的两个旧版本
  * v1.1.0 = tag a30b099（GitHub 已发布版本，origin/main 跟踪点）
  * 8f66c42 = 本轮工作之前 GitHub 主分支状态（本地 commit 1 的父提交）
  两者用 git archive 导出干净源码树后 MSBuild Release x64 完整重编译（退出码 0）：
    mopac_cpp_v110.exe  /  mopac_cpp_old_8f66c42.exe
  git diff --name-only a30b099..8f66c42 不含任何 .cpp/.h/.vcxproj 文件
  （仅 README/LICENSE/docs/.gitignore），两版本代码相同。

结论
  两个旧版本均复现 MOZYME 优化问题（首轮梯度 721 vs 官方 3.8），并且
  MOZYME+COSMO 优化直接崩溃（exit 3221225477，输出为空）。当前版本已修复
  崩溃，但梯度错误仍在（历史遗留，非本轮引入）。
  DCART 导数在 v1.1.0、8f66c42、当前版本三者逐位相同：
  O (266.532342, -385.552987, -252.046453)，官方为 (0.050406, 0.000340, -3.101008)。

运行
  python old_probe.py <old mopac.exe> [case ...]

目录
  old_probe.py            旧版本对比脚本
  old_probe_8f66c42.log   8f66c42 运行记录
  old_probe_v110.log      v1.1.0 运行记录
  results/old_8f66c42/    8f66c42 的 MOZYME 优化输出
  results_v110/           v1.1.0 的 MOZYME 优化与 DCART 输出
  results/chk_dcart/      旧/当前/官方 DCART 导数对比
  results/chk_cosmo/      MOZYME+COSMO 优化（旧版崩溃 vs 当前/官方）
