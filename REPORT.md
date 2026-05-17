# INC-Ricci 组合核心 proof-of-concept 报告

## 1. 目标与范围

本程序实现了一个 C++17 单文件 proof-of-concept，用于测试固定三角剖分环面上的组合 Ricci-flow 风格半径优化核心。

程序关注的是 INC-Ricci 思路中的快速组合部分：

- 生成周期三角剖分环面；
- 通过 2-2 edge flip 随机化拓扑；
- 为每个顶点分配圆半径 `r_i = exp(u_i)`；
- 用相切圆三角形边长 `l_ij = r_i + r_j` 计算每个三角形内角；
- 计算离散曲率 `K_i = 2*pi - sum(theta_i)`；
- 通过 Ricci-flow 更新或有限差分梯度下降使 `K_i -> 0`；
- 输出曲率、半径、拓扑和度分布诊断。

本程序不实现 Nature Materials SMC + FIRE 工作流，也不实现 PRL 级别的坐标重构 circle packing solver。

## 2. 编译与运行

编译命令：

```bash
g++ -O3 -std=c++17 inc_ricci.cpp -o inc_ricci
```

测试命令：

```bash
./inc_ricci --test
```

正式运行命令：

```bash
./inc_ricci --nx 16 --ny 16 --flips 2000 --seed 1 --max_iter 200000 --tol 1e-10 --method ricci --out run_16_ricci
./inc_ricci --nx 12 --ny 12 --flips 1000 --seed 2 --max_iter 20000 --tol 1e-8 --method grad --out run_12_grad
```

## 3. 实现概况

源码文件为 `inc_ricci.cpp`，使用 C++17 标准库，无外部依赖。

主要数据结构：

- `Face { int a,b,c; }`：三角面；
- `EdgeKey { int u,v; }`：归一化无向边，保证 `u < v`；
- `State`：保存当前半径、曲率、面角、能量与诊断量；
- `RunResult`：保存优化结果、统计历史和停止原因。

主要函数：

- `generate_periodic_triangular_torus()`：生成周期三角剖分环面；
- `build_edges()`：构建 edge-to-face 映射；
- `validate_triangulation()`：验证闭环面三角剖分；
- `attempt_random_flip()`：尝试一次合法 2-2 翻边；
- `compute_angles_for_face()`：计算相切圆三角形三个内角；
- `compute_state()` / `compute_curvature()`：计算曲率和能量；
- `normalize_u()`：保持 `sum_i u_i = 0`；
- `run_ricci()`：组合 Ricci-flow 更新；
- `run_grad()`：有限差分梯度下降调试优化器；
- `write_outputs()`：写出 CSV、faces 和 summary。

## 4. 数学模型

每个顶点有正半径：

```text
r_i = exp(u_i)
```

每条边的圆填充边长为：

```text
l_ij = r_i + r_j
```

对三角形 `(i,j,k)`，顶点 `i` 的角为：

```text
theta_i = acos(((r_i+r_j)^2 + (r_i+r_k)^2 - (r_j+r_k)^2)
               / (2*(r_i+r_j)*(r_i+r_k)))
```

离散曲率：

```text
K_i = 2*pi - sum_incident_faces theta_i
```

环面目标平均曲率为零，因此优化目标是：

```text
K_i -> 0 for all i
E_K = 0.5 * sum_i K_i^2 -> 0
```

## 5. 拓扑验证

程序在初始生成和翻边后验证：

- 每个面有三个不同顶点；
- 每条边恰好属于两个面；
- Euler characteristic `V - E + F = 0`；
- 闭三角剖分环面满足 `F = 2V` 和 `E = 3V`。

正式运行中：

| Run | V | E | F | Euler characteristic | flips requested | flips accepted |
|---|---:|---:|---:|---:|---:|---:|
| `run_16_ricci` | 256 | 768 | 512 | 0 | 2000 | 1618 |
| `run_12_grad` | 144 | 432 | 288 | 0 | 1000 | 807 |

## 6. 自动测试结果

`./inc_ricci --test` 结果：

| Test | 内容 | 结果 | 关键数值 |
|---|---|---|---|
| Test 1 | `8x8` 规则三角环面，`u_i=0` | PASS | `max_abs_K = 1.776357e-15` |
| Test 2 | `8x8`，200 次翻边，Ricci flow | PASS | `final_max_abs_K = 8.686049e-09` |
| Test 3 | `8x8`，200 次翻边，grad | PASS | `E_K` 单调下降，`final_E_K = 3.514299e-05` |

Test 1 说明规则三角环面初始态满足六个等边三角形围绕每个顶点，曲率在机器精度附近为零。

## 7. 正式运行结果

### 7.1 `run_16_ricci`

命令：

```bash
./inc_ricci --nx 16 --ny 16 --flips 2000 --seed 1 --max_iter 200000 --tol 1e-10 --method ricci --out run_16_ricci
```

结果：

| 指标 | 数值 |
|---|---:|
| iterations | 563 |
| initial `E_K` | `1.53583416110313328e+03` |
| initial `max_abs_K` | `1.80000261249385396e+01` |
| final `E_K` | `1.55527529116044540e-19` |
| final `max_abs_K` | `9.80160397290319452e-11` |
| final `r_min` | `4.02476805295839679e-02` |
| final `r_max` | `1.28835142026530054e+01` |
| `E_K` monotonic | yes |
| stop reason | `converged: max_abs_K < tol` |

结论：Ricci-flow 方法在该随机化 `16x16` 三角剖分环面上收敛到目标容差以下，得到组合意义上的零曲率相切圆三角形度量。

### 7.2 `run_12_grad`

命令：

```bash
./inc_ricci --nx 12 --ny 12 --flips 1000 --seed 2 --max_iter 20000 --tol 1e-8 --method grad --out run_12_grad
```

结果：

| 指标 | 数值 |
|---|---:|
| iterations | 1104 |
| initial `E_K` | `8.56344221163784255e+02` |
| initial `max_abs_K` | `2.29750754527028462e+01` |
| final `E_K` | `1.99067653151959583e+01` |
| final `max_abs_K` | `6.28318519052693247e+00` |
| final `r_min` | `3.12178994377280278e-02` |
| final `r_max` | `5.68950248610782504e+08` |
| `E_K` monotonic | yes |
| stop reason | `stalled: alpha < 1e-14` |

结论：有限差分梯度法保持了能量单调下降，但该运行没有收敛到零曲率。最终存在约 `2*pi` 的最大曲率残差，并伴随极端半径比，说明 steepest descent 进入了半径退化边界附近。

## 8. 额外诊断：同一拓扑 Ricci 探针

为判断 `run_12_grad` 的失败是否来自拓扑本身，额外使用同一参数但改用 Ricci-flow：

```bash
./inc_ricci --nx 12 --ny 12 --flips 1000 --seed 2 --max_iter 50000 --tol 1e-10 --method ricci --out run_12_ricci_probe
```

结果：

| 指标 | 数值 |
|---|---:|
| final `E_K` | `7.50710141000150704e-20` |
| final `max_abs_K` | `9.63980006929432420e-11` |
| final `r_min` | `3.56351200324417777e-02` |
| final `r_max` | `1.33380142741710142e+01` |
| `E_K` monotonic | yes |
| stop reason | `converged: max_abs_K < tol` |

这说明 `run_12_grad` 的未收敛不是因为随机拓扑没有零曲率解，而是有限差分梯度下降在该问题上数值路径较差，容易走向半径极端化。

## 9. 强 INC residual 解释

在本组合模型中，每个三角面都被定义为互相相切的圆三角形，边长固定为：

```text
l_ij = r_i + r_j
```

因此局部理想 steric angle 与实际三角形角度由同一个相切三角形构造给出：

```text
theta_residual = 0 by construction
```

真正需要优化的是全局闭合一致性，即每个顶点的角和是否等于 `2*pi`：

```text
curvature_residual = max_abs_K
```

运行结果：

| Run | `theta_residual` | `curvature_residual` |
|---|---:|---:|
| `run_16_ricci` | 0 by construction | `9.80160397290319452e-11` |
| `run_12_grad` | 0 by construction | `6.28318519052693247e+00` |
| `run_12_ricci_probe` | 0 by construction | `9.63980006929432420e-11` |

解释：

- `theta_residual = 0` 表示每个局部三角形都是相切圆几何，局部 steric 条件完美满足；
- `curvature_residual -> 0` 表示这些局部相切三角形可以在组合意义上闭合成平坦环面 circle-packing metric；
- 两者同时成立时，可视为本组合模型中的 strong-INC endpoint。

## 10. 输出文件

每次运行输出目录包含：

- `stats.csv`：迭代历史，包含 `iter,E_K,max_abs_K,step,r_min,r_max,accepted`；
- `radii.csv`：最终 `u_i` 与 `r_i`；
- `faces.dat`：最终三角面顶点编号；
- `degree_hist.csv`：顶点度分布；
- `summary.txt`：参数、拓扑验证、最终诊断、停止原因和限制说明。

## 11. 结论

该 C++17 proof-of-concept 已实现并验证了 INC-Ricci 的组合半径/曲率核心。

最重要的结果是：在随机翻边后的固定三角剖分环面上，Ricci-flow 更新能够快速、单调地把离散曲率残差压到 `1e-10` 量级，说明组合相切圆三角形可以被驱动到平坦环面度量。

有限差分梯度法作为调试优化器可以保证能量单调下降，但对较困难随机拓扑可能走向半径退化边界，不应作为主要生成器。实际 PoC 中应优先使用 `--method ricci`。

## 12. 限制

本程序只测试 INC-Ricci 思路中的组合 radius/curvature 部分：

- 不重构圆心的显式欧氏坐标；
- 不验证坐标级 circle packing 嵌入；
- 不实现 Nature Materials SMC + FIRE protocol；
- 不证明物理 protocol 的收敛性；
- 只验证随机化三角剖分环面是否能被驱动到零曲率相切圆组合度量。
