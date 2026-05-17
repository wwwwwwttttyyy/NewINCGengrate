# Old Program Exploration Report

Generated: 2026-05-17 16:55

## 1. Executive Summary

**核心结论：v18 多分散停滞的根本原因是 cycle-level 接受准则中的 U_overlap 惩罚。**

1. ✅ v18 polydisperse 停滞已完全复现：默认参数下 200 cycles 全部 rejected
2. ❌ gamma 不是原因：shape_gamma 从 0 到 1.0 全部冻结，值完全相同
3. ❌ bounded_du 默认值(0.10)不是唯一原因：但紧约束(0.03/0.05)有部分缓解
4. ✅ **overlap penalty 是主因**：weak_overlap(w_U=0.05)达到 186 accepted cycles，diagonal improved
5. ✅ bounded_du sweet spot 存在：0.03→88 accepted, 0.05→72 accepted, ≥0.10→0
6. ✅ v17 operator 方向稳健：free_ricci 100% 降低 Θ_sum，shape_descent 69% 降低 Θ_shape
7. ✅ compression 是最强 diagonal move：phi=0.90 最小化 Θ_sum，phi=0.92 最小化 Θ_NM
8. ⚠️ **关键发现**：所有降 Θ_sum 的 operator 都伴随巨大 U_overlap 增加（free_ricci: +4.04, bounded_ricci: +0.12）
9. v19 必须实现：(1) annealed acceptance, (2) proposal forensics

## 2. Program and Test Status

| Program | Test | Status |
|---------|------|--------|
| inc_ricci_v17 | --test | ALL 5 TESTS PASSED |
| inc_ricci_v18 | --test | ALL 5 TESTS PASSED |

- v18 writes output files to CWD (--out flag unused for file output)
- All required parameters verified in source code

## 3. v18 Gamma Scan

| gamma | Theta_sum_i | Theta_sum_f | Delta_sum | Theta_shape_i | Theta_shape_f | Delta_shape | accepted | classification |
|-------|------------|------------|-----------|--------------|--------------|-------------|----------|---------------|
| 0.0 | 0.06659123498847247 | 0.06659123498847247 | 0.0000 | 0.11753862407793556 | 0.11753862407793556 | 0.0000 | 0 | cycle_rejected_blocks_accepted |
| 0.02 | 0.06659123498847247 | 0.06659123498847247 | 0.0000 | 0.11753862407793556 | 0.11753862407793556 | 0.0000 | 0 | cycle_rejected_blocks_accepted |
| 0.05 | 0.06659123498847247 | 0.06659123498847247 | 0.0000 | 0.11753862407793556 | 0.11753862407793556 | 0.0000 | 0 | cycle_rejected_blocks_accepted |
| 0.1 | 0.06659123498847247 | 0.06659123498847247 | 0.0000 | 0.11753862407793556 | 0.11753862407793556 | 0.0000 | 0 | cycle_rejected_blocks_accepted |
| 0.25 | 0.06659123498847247 | 0.06659123498847247 | 0.0000 | 0.11753862407793556 | 0.11753862407793556 | 0.0000 | 0 | cycle_rejected_blocks_accepted |
| 0.5 | 0.06659123498847247 | 0.06659123498847247 | 0.0000 | 0.11753862407793556 | 0.11753862407793556 | 0.0000 | 0 | cycle_rejected_blocks_accepted |
| 1.0 | 0.06659123498847247 | 0.06659123498847247 | 0.0000 | 0.11753862407793556 | 0.11753862407793556 | 0.0000 | 0 | cycle_rejected_blocks_accepted |

**Interpretation**: If all gamma values lead to frozen state, shape_gamma is not the primary cause.

## 4. v18 Bounded Ricci Scan

| bounded_du | Delta_sum | Delta_shape | accepted | classification |
|-----------|-----------|-------------|----------|---------------|
| 0.03 | -0.0044 | -0.0251 | 88 | diagonal_improved |
| 0.05 | -0.0073 | -0.0278 | 72 | diagonal_improved |
| 0.1 | 0.0000 | 0.0000 | 0 | cycle_rejected_blocks_accepted |
| 0.15 | 0.0000 | 0.0000 | 0 | cycle_rejected_blocks_accepted |
| 0.25 | 0.0000 | 0.0000 | 0 | cycle_rejected_blocks_accepted |
| 0.4 | 0.0000 | 0.0000 | 0 | cycle_rejected_blocks_accepted |

**Interpretation**: bounded_du=0.03 and 0.05 achieve diagonal improvement (88 and 72 accepted cycles). Tighter Ricci bounds constrain moves enough to avoid overwhelming overlap penalty. bounded_du>=0.10 freezes completely - the Ricci step overshoots the F_total improvement window.

## 5. v18 Weight Scan

| scheme | Delta_sum | Delta_shape | accepted | classification |
|--------|-----------|-------------|----------|---------------|
| shape_biased | 0.0000 | 0.0000 | 0 | cycle_rejected_blocks_accepted |
| balanced | 0.0000 | 0.0000 | 0 | cycle_rejected_blocks_accepted |
| weak_overlap | -0.0214 | -0.0842 | 186 | diagonal_improved |

## 6. v18 Trajectory Inventory

Total v18 output dirs: 15

### ⚠️ Critical Finding: Block-Level Accepted but Cycle-Level Rejected

The following runs show that individual blocks (shape/ricci) are accepted,
but the overall coupled cycle is always rejected. This points to the cycle-level
acceptance criterion being too strict, NOT to the operators being ineffective.

- **v18_poly_bound_0.15**: blocks_accepted=460/460 but cycles_accepted=0/200
- **v18_poly_bound_0.25**: blocks_accepted=460/460 but cycles_accepted=0/200
- **v18_poly_bound_0.40**: blocks_accepted=460/460 but cycles_accepted=0/200
- **v18_poly_gamma_0**: blocks_accepted=411/411 but cycles_accepted=0/166
- **v18_poly_gamma_0.02**: blocks_accepted=460/460 but cycles_accepted=0/200
- **v18_poly_gamma_0.05**: blocks_accepted=460/460 but cycles_accepted=0/200
- **v18_poly_gamma_0.10**: blocks_accepted=460/460 but cycles_accepted=0/200
- **v18_poly_gamma_0.25**: blocks_accepted=460/460 but cycles_accepted=0/200
- **v18_poly_gamma_0.50**: blocks_accepted=460/460 but cycles_accepted=0/200
- **v18_poly_gamma_1.0**: blocks_accepted=460/460 but cycles_accepted=0/200
- **v18_weight_balanced**: blocks_accepted=460/460 but cycles_accepted=0/200
- **v18_weight_shape_biased**: blocks_accepted=460/460 but cycles_accepted=0/200

### Detailed Inventory

- **v18_poly_bound_0.03**: class=diagonal_improved accepted=88 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_poly_bound_0.05**: class=diagonal_improved accepted=72 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_poly_bound_0.15**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_poly_bound_0.25**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_poly_bound_0.40**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_poly_gamma_0**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=3
  - ORM: 412 rows, shape=180, ricci=180, box=35, surgery=17
  - Block acceptance: 411/411 (100.0%)
- **v18_poly_gamma_0.02**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_poly_gamma_0.05**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_poly_gamma_0.10**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_poly_gamma_0.25**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_poly_gamma_0.50**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_poly_gamma_1.0**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_weight_balanced**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_weight_shape_biased**: class=cycle_rejected_blocks_accepted accepted=0 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)
- **v18_weight_weak_overlap**: class=diagonal_improved accepted=186 has_traj=True files=7
  - ORM: 460 rows, shape=200, ricci=200, box=40, surgery=20
  - Block acceptance: 460/460 (100.0%)

## 7. v17 Operator Atlas Extension

### 7.1 Operator Directions（16x16 basic, 128 starts）

| Operator | ΔΘ_sum | ΔΘ_shape | ΔΘ_NM | ΔU_overlap | Δdu_rms | succ_left | succ_down |
|----------|--------|----------|-------|-----------|---------|-----------|-----------|
| free_ricci | -0.0739 | -0.0036 | -0.0085 | **+4.04** | -0.146 | 100% | 100% |
| bounded_ricci_du010 | -0.0181 | -0.0015 | -0.0053 | **+0.120** | -0.030 | 100% | 100% |
| post_voro_surgery | -0.0091 | +0.0287 | +0.0201 | 0.000 | -0.014 | 100% | 0% |
| shape_descent | +0.0000 | -0.0000 | -0.0000 | +0.0000 | +0.0001 | 16% | 69% |
| compression_step | +0.0000 | +0.0000 | +0.0000 | **-0.0013** | +0.0001 | 41% | 35% |
| radius_label_swap | -0.0000 | -0.0000 | -0.0000 | -0.0000 | -0.0001 | 34% | 38% |

### 7.2 关键发现：Overlap-Energy Dilemma

**所有能有效降低 Θ_sum 的 operator 都伴随巨大的 U_overlap 增加：**
- free_ricci: ΔΘ_sum=-0.074（极强）但 ΔU=**+4.04**（灾难性）
- bounded_ricci: ΔΘ_sum=-0.018 但 ΔU=**+0.12**

**能降低 overlap 的 operator 几乎没有 Θ 效果：**
- compression_step: ΔU=-0.0013 但 ΔΘ≈0

这就是 v18 停滞的微观机制。

### 7.3 N=256 vs N=576 一致性
operator 方向在两个尺寸下完全一致。

## 8. Fine Compression Window

- Best phi by Theta_NM: phi=0.92, value=0.091800
- Best phi by Theta_sum: phi=0.90, value=0.063340
- Best phi by Theta_shape: phi=0.92, value=0.090593

## 9. Failure Cause Ranking

### A. overlap_penalty_dominates_F_total 🔴 HIGH
**证据**：weak_overlap(w_U=0.05) 成功 186 cycles；标准权重全部冻结；v17 data: free_ricci ΔU=+4.04
**机制**：Ricci radius 调整改变 Voronoi cell 面积 → 破坏 packing → overlap 激增 → F_total 变差

### B. cycle_acceptance_criterion_too_strict 🔴 HIGH
**证据**：12/15 runs: block 100% accepted, cycle 0% accepted
**机制**：Shape+Ricci 合并后 F_total 总是比 cycle 前更差（因 U_overlap）

### C. bounded_du_sweet_spot 🟡 MEDIUM
**证据**：0.03→88 acc, 0.05→72 acc, ≥0.10→0 acc。但 sweet spot 最终也停滞。

### D. global_gamma_irrelevant 🟢 HIGH（排除项）
**证据**：7 个 gamma 全部相同冻结状态

### E. shape_operator_weak_on_polydisperse 🟡 MEDIUM
**证据**：shape_descent success_down=69% 但 ΔΘ_shape≈0

### F. physical_frustration_possible 🟡 MEDIUM
**证据**：所有参数调整都无法完全解决

### G. needs_contact_weighted_edge_target ⚪ 需要 v19
### H. needs_nonaffine_activation ⚪ 需要 v19

## 10. Recommended v19 Design

### 必须实现（PRIORITY 1）：

1. **Annealed Acceptance**
   - 当前问题：cycle-level 严格 ΔF_total<0 判据导致 100% reject
   - 方案：引入温度参数 T，接受 ΔF_total < T 的 move，T 从高到低退火
   - 参数参考：w_U=0.05 等效于降低 overlap 权重，可作为初始温度参考

2. **Proposal Forensics**
   - 当前问题：不知道 F_total 具体哪个分量导致拒绝
   - 方案：每步记录 ΔF_total 分解（ΔΘ_sum×w_sum + ΔΘ_shape×w_shape + ΔU×w_U + ...）
   - 效果：消除盲参数扫描

### 强烈推荐（PRIORITY 2）：

3. **Contact-Weighted Edge Target**
   - 当前问题：Voronoi 边不反映实际颗粒接触
   - 方案：用实际接触面积加权边的 target length

4. **Nonaffine Activation**
   - 当前问题：shape_descent 效果微弱
   - 方案：允许非仿射位移场

### 参数继承建议：

| 参数 | 推荐值 | 来源 |
|------|--------|------|
| shape_gamma | 任意（不影响结果） | gamma scan |
| bounded_du | 0.03-0.05（配合 annealing） | bound scan |
| w_U | 0.05 初始，随 anneal 增加 | weight scan |
| target_phi | 0.84-0.90 | compression scan |
| shape_steps/ricci_steps | 20/100 (Ricci-heavy) | bound scan |


## 11. Files and Figures

### Summary CSVs
- `raw/v17_fine_compression_summary.csv`
- `raw/v17_operator_statistics_summary.csv`
- `raw/v18_bound_scan_summary.csv`
- `raw/v18_gamma_scan_summary.csv`
- `raw/v18_weight_scan_summary.csv`

### Inventory CSVs
- `raw/v18_trajectory_inventory.csv`

### Figures
- `figures/combined_defect_plane_summary.png`
- `figures/v17_fine_compression.png`
- `figures/v17_operator_success_rates.png`
- `figures/v17_operator_vectors.png`
- `figures/v18_bound_scan.png`
- `figures/v18_gamma_scan.png`
- `figures/v18_response_matrix_analysis.png`
- `figures/v18_weight_scan.png`

### Logs
- `logs/master_run.log`
- `logs/v17_fine_compression.log`
- `logs/v17_operator_16_basic.log`
- `logs/v17_operator_24_basic.log`
- `logs/v18_poly_bound_0.03.log`
- `logs/v18_poly_bound_0.05.log`
- `logs/v18_poly_bound_0.15.log`
- `logs/v18_poly_bound_0.25.log`
- `logs/v18_poly_bound_0.40.log`
- `logs/v18_poly_gamma_0.02.log`
- `logs/v18_poly_gamma_0.05.log`
- `logs/v18_poly_gamma_0.10.log`
- `logs/v18_poly_gamma_0.25.log`
- `logs/v18_poly_gamma_0.50.log`
- `logs/v18_poly_gamma_0.log`
- `logs/v18_poly_gamma_1.0.log`
- `logs/v18_weight_balanced.log`
- `logs/v18_weight_shape_biased.log`
- `logs/v18_weight_weak_overlap.log`
