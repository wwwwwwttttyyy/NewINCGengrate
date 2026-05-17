# INC-Ricci v12 coordinate-optimized voro++ topology report

## Build And Tests

Build command used under WSL:

```bash
g++ -O3 -std=c++17 inc_ricci_v12.cpp -o inc_ricci_v12 \
  -I/home/tianyu/miniconda3/envs/analyse/include \
  -I/home/tianyu/miniconda3/envs/analyse/include/voro++ \
  -L/home/tianyu/miniconda3/envs/analyse/lib \
  -lvoro++ \
  -Wl,-rpath,/home/tianyu/miniconda3/envs/analyse/lib \
  -pthread
```

Test command:

```bash
./inc_ricci_v12 --test
```

Result: all v12 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Triangular lattice equal radii | PASS | `B_rms=1.776357e-15`, `du_rms=2.761069e-22`, `E=192`, `F=128` |
| Soft relaxation | PASS | `U_initial=4.534427`, `U_final=0.0102323` |
| Coordinate evaluation, N=64 | PASS | `S_pred=120.2802`, `du_rms=0.1713305` |
| Small coordinate optimization | PASS | `S_X=55.94976 -> 55.94672` |
| Coordinate optimization strategy compare | PASS | 6 strategy rows written |

Implementation note: v12 uses two-stage coordinate scoring. It evaluates all coordinate candidates with a cheap score based on voro++ topology, angle budget, and overlap. For performance, full finite-difference linear-response evaluation is done periodically and at the final state. The final reported `S_pred`, `du_rms`, and `du_max` are full predictor values.

## v12 Implementation Summary

`inc_ricci_v12.cpp` preserves the v11 voro++ radical topology path and adds:

- `evaluate_coordinates(X,r)` using voro++ radical cells, dual triangulation, angle-budget diagnostics, overlap energy, and optional full linear predictor;
- coordinate optimization mode `--coord_opt 1`;
- move modes: `single_particle`, `local_cluster`, `swap_positions`, `soft_relax_step`, `mixed`;
- two-stage scoring with `S_X`;
- `coord_opt_stats.csv`, `coords_optimized.csv`, `coord_opt_summary.txt`, `voro_delaunay_summary_final.txt`, and `linear_response_summary_final.txt`;
- coordinate optimization strategy comparison;
- coordinate optimization multistart;
- optional phi scan.

Score:

```text
S_X =
w_pred   * log10(1 + S_pred / S_pred_scale)
+ w_B     * log10(1 + B_rms / B_scale)
+ w_durms * (du_rms / du_rms_target)^2
+ w_dumax * (du_max / du_max_target)^2
+ w_ov    * log10(1 + U_overlap / U_scale)
+ w_inv   * invalid_topology_penalty
```

## Single N=256 Coordinate Optimization

Command:

```bash
./inc_ricci_v12 --nx 16 --ny 16 \
  --radius_seed 101 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --coord_mode soft_relax \
  --target_phi 0.84 \
  --relax_steps 30000 \
  --use_voropp 1 \
  --coord_opt 1 \
  --coord_opt_steps 3000 \
  --coord_opt_trials 20 \
  --coord_opt_move_scale 0.02 \
  --coord_move_mode mixed \
  --two_stage_eval 1 \
  --coord_full_eval_top 3 \
  --out coord_opt_v12_16_delta025_ratio5 \
  --progress_every 100
```

| Diagnostic | Initial | Final |
|---|---:|---:|
| `S_X` | `77.526811` | `77.234662` |
| `B_rms` | `0.528627` | `0.519804` |
| `B_max` | `1.896204` | `1.896204` |
| `S_pred` | `129.924895` | `129.400383` |
| `du_rms` | `0.140375` | `0.139132` |
| `du_max` | `0.615459` | `0.616352` |
| `U_overlap` | `1.6315e-04` | `1.5774e-04` |
| `max_overlap` | `2.1699e-03` | `2.1105e-03` |
| degree-radius corr | `0.790293` | `0.795733` |
| valid torus | yes | yes |
| `linear_predicted_good` | no | no |

Accepted coordinate moves: `2972 / 3000`.

Interpretation: mixed coordinate moves lowered the coordinate score slightly and improved `B_rms`, `S_pred`, and `du_rms` marginally. The main obstruction did not change: `du_max` stayed around `0.62`, well above the physical threshold `0.1`.

## Strategy Compare, N=256

Command:

```bash
./inc_ricci_v12 --nx 16 --ny 16 \
  --seed 1 \
  --radius_seed 101 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --coord_mode soft_relax \
  --target_phi 0.84 \
  --relax_steps 30000 \
  --use_voropp 1 \
  --coord_opt_strategy_compare 1 \
  --coord_opt_steps 2000 \
  --coord_opt_trials 20 \
  --two_stage_eval 1 \
  --coord_full_eval_top 3 \
  --out coord_opt_strategy_v12_16 \
  --progress_every 0
```

| Strategy | `S_X` | `B_rms` | `S_pred` | `du_rms` | `du_max` | max overlap | degree corr | predicted |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `soft_relax_only` | `77.5268` | `0.5286` | `129.9249` | `0.1404` | `0.6155` | `0.00217` | `0.7903` | no |
| `coord_opt_single_particle` | `61.0752` | `0.4974` | `104.0486` | `0.1325` | `0.5228` | `0.00476` | `0.8070` | no |
| `coord_opt_local_cluster` | `77.5268` | `0.5286` | `129.9249` | `0.1404` | `0.6155` | `0.00217` | `0.7903` | no |
| `coord_opt_swap` | `77.5268` | `0.5286` | `129.9249` | `0.1404` | `0.6155` | `0.00217` | `0.7903` | no |
| `coord_opt_mixed` | `77.2396` | `0.5198` | `129.4004` | `0.1391` | `0.6164` | `0.00213` | `0.7957` | no |
| `coord_opt_mixed_plus_post_surgery` | `33.3254` | `0.4660` | `61.9257` | `0.1220` | `0.3037` | `0.00213` | `0.8116` | no |

Best pure coordinate move: `single_particle`.

Best overall: `coord_opt_mixed_plus_post_surgery`, but this includes post-coordinate topology surgery, so the topology is no longer exact radical Delaunay. It improves `du_rms` and `du_max` substantially, but still misses the strict thresholds.

## N=256 Coordinate-Optimization Multistart

Command:

```bash
./inc_ricci_v12 --nx 16 --ny 16 \
  --seed 1000 \
  --radius_seed 101 \
  --same_radii_all_starts 1 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --coord_mode soft_relax \
  --target_phi 0.84 \
  --relax_steps 20000 \
  --use_voropp 1 \
  --coord_opt 1 \
  --coord_opt_steps 2000 \
  --coord_opt_trials 16 \
  --coord_move_mode mixed \
  --two_stage_eval 1 \
  --coord_full_eval_top 2 \
  --coord_opt_multistart 1 \
  --threads 16 \
  --multistart_count 32 \
  --coord_opt_multistart_out coord_opt_multistart_v12_16.csv \
  --progress_every 0
```

| Statistic | Value |
|---|---:|
| rows | 32 |
| valid torus rate | `32 / 32` |
| `linear_predicted_good` | `0 / 32` |
| best `du_rms` | `0.1298103177` |
| median `du_rms` | `0.1452826070` |
| best `du_max` | `0.3918705741` |
| median `du_max` | `0.5125044630` |
| best `S_pred` | `81.95102676` |
| median `S_pred` | `108.9316844` |
| starts beating v11 best `du_rms ~= 0.1298103` | `1 / 32` |

Best by `du_rms`:

| Field | Value |
|---|---:|
| run id | `25` |
| `S_pred_final` | `100.7850568` |
| `du_rms_final` | `0.1298103177` |
| `du_max_final` | `0.5146946675` |
| degree corr | `0.8216611406` |

Best by `S_pred`:

| Field | Value |
|---|---:|
| run id | `9` |
| `S_pred_final` | `81.95102676` |
| `du_rms_final` | `0.1346863008` |
| `du_max_final` | `0.3942882449` |

Interpretation: coordinate optimization did not materially beat the best v11 multistart result in `du_rms`, but it found lower `du_max` in the best-score run.

## N=576 Coordinate-Optimization Multistart

Command:

```bash
./inc_ricci_v12 --nx 24 --ny 24 \
  --seed 2000 \
  --radius_seed 202 \
  --same_radii_all_starts 1 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --coord_mode soft_relax \
  --target_phi 0.84 \
  --relax_steps 30000 \
  --use_voropp 1 \
  --coord_opt 1 \
  --coord_opt_steps 1500 \
  --coord_opt_trials 12 \
  --coord_move_mode mixed \
  --two_stage_eval 1 \
  --coord_full_eval_top 2 \
  --coord_opt_multistart 1 \
  --threads 16 \
  --multistart_count 12 \
  --coord_opt_multistart_out coord_opt_multistart_v12_24.csv \
  --progress_every 0
```

| Statistic | Value |
|---|---:|
| rows | 12 |
| valid torus rate | `12 / 12` |
| `linear_predicted_good` | `0 / 12` |
| best `du_rms` | `0.1384877873` |
| median `du_rms` | `0.1449247461` |
| best `du_max` | `0.4504609373` |
| median `du_max` | `0.5317222682` |
| best `S_pred` | `97.17293218` |
| median `S_pred` | `114.8707812` |
| runtime | about `519 s` |

The larger system shows the same qualitative behavior: voro++ topology stays valid, but ordinary coordinate moves do not reach small-du closure.

## Phi Scan, N=256

Command used `phi_list=0.76,0.78,0.80,0.82,0.84,0.86`.

| phi | `S_pred` | `du_rms` | `du_max` | `B_rms` | max overlap | predicted |
|---:|---:|---:|---:|---:|---:|---:|
| 0.76 | `204.1110` | `0.2009` | `0.7201` | `0.7139` | `3.85e-08` | no |
| 0.78 | `183.1375` | `0.2073` | `0.6206` | `0.7731` | `5.15e-08` | no |
| 0.80 | `144.0707` | `0.1795` | `0.5567` | `0.6189` | `1.71e-07` | no |
| 0.82 | `138.3384` | `0.1637` | `0.5842` | `0.6072` | `8.39e-08` | no |
| 0.84 | `106.4325` | `0.1444` | `0.4980` | `0.5380` | `0.00388` | no |
| 0.86 | `102.0492` | `0.1398` | `0.4912` | `0.5083` | `0.05846` | no |

Higher density improved the linear predictor but increased residual overlaps. No phi produced `linear_predicted_good`.

## Conclusion

v12 tests the next hypothesis directly: can fixed-radius real-space coordinate moves improve the induced radical topology enough to reach small-radius-perturbation closure?

Result:

```text
coordinate optimization improves S_X: slightly
single-particle coordinate moves help: yes
mixed coordinate moves help: weakly
post-coordinate topology surgery helps: strongly, but topology is no longer exact radical Delaunay
strict physical threshold reached: no
```

The best pure coordinate optimization remained around:

```text
du_rms ~= 0.13
du_max ~= 0.39-0.52
```

This is better than random coordinate topology, but still far from:

```text
du_rms < 0.03
du_max < 0.1
```

Interpretation:

- The v11/v12 voro++ topology generator is robust; topology validity is no longer the bottleneck.
- Ordinary small coordinate moves mostly smooth overlap and sometimes change topology beneficially, but they do not reliably enter a new small-du basin.
- Post-voro or post-coordinate topology surgery remains the strongest improvement mechanism, suggesting that coordinate optimization and controlled topology surgery are complementary.
- If the physical INC-like regime exists at `Delta=0.25`, `radius_ratio=5`, the missing process is probably not simple overlap relaxation. It is more likely a stronger coordinate/topology equilibration process, explicit swap/annealing, or direct optimization of angle-budget/linear-response objectives through larger topology-changing coordinate moves.
