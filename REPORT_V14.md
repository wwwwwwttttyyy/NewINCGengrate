# INC-Ricci v14 soft-hard diagnostic report

## Build And Tests

Build command used under WSL:

```bash
g++ -O3 -std=c++17 inc_ricci_v14.cpp -o inc_ricci_v14 \
  -I/home/tianyu/miniconda3/envs/analyse/include \
  -I/home/tianyu/miniconda3/envs/analyse/include/voro++ \
  -L/home/tianyu/miniconda3/envs/analyse/lib \
  -lvoro++ \
  -Wl,-rpath,/home/tianyu/miniconda3/envs/analyse/lib \
  -pthread
```

Test command:

```bash
./inc_ricci_v14 --test
```

Result: all v14 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Rectangular hex equal-radii baseline | PASS | `Theta_ring=3.47e-16`, `B_rms=1.78e-15`, `du_rms=2.76e-22`, `E=768`, `F=512` |
| Face/ring theta comparison on hex baseline | PASS | `Theta_face=2.13e-16`, `Theta_ring=2.13e-16`, difference `0` |
| Soft polydisperse N=64 diagnostics | PASS | `Theta=0.129382`, `B_rms=0.543094`, `du_rms=0.137687` |
| Controlled perturb suite small | PASS | `Theta` increased from `2.13e-16` to `3.19e-02` |
| Radius escape small | PASS | stats file produced, final ratio `3.447` |

## v14 Implementation Summary

`inc_ricci_v14.cpp` preserves the v13/v12 voro++ radical topology and linear-response pipeline, and adds:

- rectangular hexagonal box mode fixing the equal-radii `Theta_NM` zero baseline;
- face-based and neighbor-ring `Theta_NM` definitions;
- primary `Theta_NM = Theta_ring_mean_abs`;
- soft-hard classification: `soft_low`, `hard_low`, and four regimes;
- controlled perturbation suite;
- soft-hard multistart scan;
- sample-class comparison;
- radius-escape diagnostic;
- scan-level soft-hard correlation summaries.

## Baseline Hex

Command:

```bash
./inc_ricci_v14 --nx 16 --ny 16 --box_mode rectangular_hex --baseline_hex_test 1 --compute_nm_theta 1 --out baseline_hex_v14_16 --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| `Theta_face_mean_abs` | `3.469446951953614e-16` |
| `Theta_ring_mean_abs` | `3.469446951953614e-16` |
| `B_rms` | `1.776356839400250e-15` |
| `du_rms` | `2.761582819331815e-22` |
| `du_max` | `2.851910971030779e-22` |

The zero baseline is fixed. Absolute `Theta_NM` values are now interpretable.

## N=256 Soft-Hard Diagnostic

Command:

```bash
./inc_ricci_v14 --nx 16 --ny 16 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --box_mode square --target_phi 0.84 --relax_steps 30000 --use_voropp 1 --compute_nm_theta 1 --out soft_hard_diag_v14_16_delta025_ratio5 --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| `Theta_NM` | `0.1181521953` |
| `Theta_face_ring_difference` | `2.498e-16` |
| `B_rms` | `0.5286270607` |
| `S_pred` | `129.9248952` |
| `rms_K_lin` | `8.164975892e-04` |
| `du_rms` | `0.1403753400` |
| `du_max` | `0.6154590162` |
| regime | `soft_high_hard_high` |

Vertex correlations:

| Correlation | Value |
|---|---:|
| `corr(Theta_i, abs(B_i))` | `0.5730815442` |
| `corr(Theta_i, abs(du_i))` | `0.5326892507` |
| `corr(Theta_i, abs(K_lin_i))` | `0.2884886315` |

## Controlled Perturbations

Hex baseline with coordinate noise:

| amp | `Theta_NM` | `B_rms` | `du_rms` | regime |
|---:|---:|---:|---:|---|
| 0 | `3.47e-16` | `1.78e-15` | `2.76e-22` | `soft_low_hard_low` |
| 0.05 | `0.03565` | `1.78e-15` | `2.76e-22` | `soft_low_hard_low` |
| 0.10 | `0.07015` | `1.78e-15` | `2.76e-22` | `soft_high_hard_low` |
| 0.20 | `0.13889` | `1.78e-15` | `2.76e-22` | `soft_high_hard_low` |

Coordinate noise increases the true coordinate-level steric defect while the radical topology remains unchanged, so the hard topology metric stays zero. This is a genuine split regime.

Hex baseline with radius noise:

| amp | `Theta_NM` | `B_rms` | `du_rms` | regime |
|---:|---:|---:|---:|---|
| 0 | `3.47e-16` | `1.78e-15` | `2.76e-22` | `soft_low_hard_low` |
| 0.02 | `0.01191` | `0.07844` | `0.02013` | `soft_low_hard_low` |
| 0.05 | `0.03069` | `0.20281` | `0.05376` | `soft_low_hard_high` |
| 0.10 | `0.06194` | `0.41823` | `0.10804` | `soft_high_hard_high` |
| 0.20 | `0.10295` | `0.68182` | `0.19091` | `soft_high_hard_high` |

For radius perturbations, soft and hard metrics scale together strongly.

## Soft-Hard Multistart, N=256

Command:

```bash
./inc_ricci_v14 --nx 16 --ny 16 --seed 1000 --radius_seed 101 --same_radii_all_starts 1 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 20000 --use_voropp 1 --compute_nm_theta 1 --soft_hard_multistart 1 --threads 16 --multistart_count 32 --soft_hard_multistart_out soft_hard_multistart_v14_16.csv --progress_every 0
```

| Statistic | Value |
|---|---:|
| rows | 32 |
| best `Theta_NM` | `0.1105432785` |
| median `Theta_NM` | `0.1197853312` |
| best `du_rms` | `0.1312777890` |
| median `du_rms` | `0.1452445940` |
| best `du_max` | `0.3918705741` |
| median `du_max` | `0.5108973635` |
| regimes | `soft_high_hard_high: 32 / 32` |

Configuration correlations:

| Correlation | Value |
|---|---:|
| `corr(Theta_NM, B_rms)` | `0.5406` |
| `corr(Theta_NM, S_pred)` | `0.3339` |
| `corr(Theta_NM, du_rms)` | `0.6093` |
| `corr(Theta_NM, du_max)` | `0.1570` |
| `corr(B_rms, du_rms)` | `0.8705` |

Interpretation: in soft-relaxed polydisperse voro samples, `Theta_NM` is a moderate proxy for `du_rms`, but not for `du_max`.

## Radius Escape

Command:

```bash
./inc_ricci_v14 --nx 16 --ny 16 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 20000 --use_voropp 1 --radius_escape_test 1 --escape_ratio_stop 1e4 --escape_max_iter 200000 --out radius_escape_v14_16_delta025_ratio5 --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| final `max_abs_K` | `8.695866e-11` |
| final radius ratio | `5.3355167247` |
| topology-distribution mismatch | no |

For this voro topology, free Ricci closure does not require extreme radius escape. It exceeds the target ratio slightly, but not catastrophically.

## Sample Class Compare

| Class | `Theta_NM` | `B_rms` | `S_pred` | `du_rms` | `du_max` | Regime |
|---|---:|---:|---:|---:|---:|---|
| `hex_equal_baseline` | `3.47e-16` | `1.78e-15` | `8.49e-13` | `2.76e-22` | `2.85e-22` | `soft_low_hard_low` |
| `voro_soft` | `0.117236` | `0.488451` | `103.608` | `0.131671` | `0.522692` | `soft_high_hard_high` |
| `abstract_v9_like` | n/a | `0.650281` | `138.414` | `0.183400` | `0.515275` | `unknown` |

The coordinate-generated voro soft sample beats the abstract sample on the hard proxy, but is not soft-low or hard-low.

## Conclusion

v14 separates the two notions cleanly:

```text
zero Theta_NM baseline fixed: yes
ordinary voro soft samples soft-low: no
ordinary voro soft samples hard-low: no
soft/hard split exists under controlled perturbations: yes
soft and hard correlated in polydisperse soft samples: moderately
```

Interpretation:

- `Theta_NM` and Ricci/linear response are related diagnostics, not identical ones.
- Coordinate perturbations of a perfect hex topology create `soft_high_hard_low`: coordinate steric order degrades while topology-level hard closure remains exact.
- Radius perturbations create strong agreement between soft and hard metrics.
- The physically relevant polydisperse voro samples remain `soft_high_hard_high`, so they are not close to either NM-like soft order or PRL/Ricci hard closure.
- Radius escape is not catastrophic for the tested voro topology, so the remaining obstruction is not simply free-Ricci singularity.
