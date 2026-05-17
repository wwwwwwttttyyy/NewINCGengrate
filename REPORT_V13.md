# INC-Ricci v13 NM-theta diagnostic report

## Build And Tests

Build command used under WSL:

```bash
g++ -O3 -std=c++17 inc_ricci_v13.cpp -o inc_ricci_v13 \
  -I/home/tianyu/miniconda3/envs/analyse/include \
  -I/home/tianyu/miniconda3/envs/analyse/include/voro++ \
  -L/home/tianyu/miniconda3/envs/analyse/lib \
  -lvoro++ \
  -Wl,-rpath,/home/tianyu/miniconda3/envs/analyse/lib \
  -pthread
```

Test command:

```bash
./inc_ricci_v13 --test
```

Result: all v13 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Triangular lattice equal radii | PASS | `B_rms=1.776357e-15`, `du_rms=2.761069e-22`, `Theta=0.0799349` |
| Soft random NM theta | PASS | `Theta=0.1267024`, `U_final=0.0102323` |
| Correlation summary | PASS | `corr(Theta,|B|)=0.445954`, `corr(Theta,|du|)=0.533859` |
| Small theta coordinate optimization | PASS | `Theta=0.145583 -> 0.145546` |
| Small theta multistart | PASS | 4 rows written |

Note: the equal-radii triangular-lattice test gives exact Ricci/linear closure, but `Theta_NM` is not machine zero because the current voro++ slab uses a square periodic box. This exposes the intended distinction between topology-level closure and coordinate-level steric geometry.

## v13 Implementation Summary

`inc_ricci_v13.cpp` preserves v12 and adds the true coordinate-level NM-style steric angle diagnostic:

- `Theta_NM` from real coordinates and lifted voro radical Delaunay faces;
- per-corner actual angle vs tangent-circle ideal angle;
- per-vertex `Theta_i`;
- vertex correlations with `|B_i|`, `|du_i|`, degree, `u0_i`, and `|K_lin_i|`;
- separate `nm_like_candidate` and `prl_like_strong_candidate` flags;
- diagnostic score `S_theta`;
- optional theta-aware coordinate optimization;
- theta multistart and theta phi scan.

Outputs added:

- `nm_theta_summary.txt`
- `nm_theta_vertex.csv`
- `nm_theta_face_corners.csv`
- `defect_correlation_summary.txt`
- `theta_coord_opt_stats.csv`
- `theta_multistart.csv`
- `theta_phi_scan.csv`

## Single N=256 Theta Diagnostic

Command:

```bash
./inc_ricci_v13 --nx 16 --ny 16 \
  --radius_seed 101 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --coord_mode soft_relax \
  --target_phi 0.84 \
  --relax_steps 30000 \
  --use_voropp 1 \
  --compute_nm_theta 1 \
  --out theta_diag_v13_16_delta025_ratio5_phi084 \
  --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| `Theta_mean_abs` | `0.1181521953` |
| `Theta_rms` | `0.1803702531` |
| `Theta_max_abs` | `0.7967878267` |
| `B_rms` | `0.5286270607` |
| `B_max` | `1.8962043324` |
| `S_pred` | `129.9248952` |
| `du_rms` | `0.1403753400` |
| `du_max` | `0.6154590162` |
| `S_theta` | `2.7419744983` |
| `nm_theta_low` | no |
| `ricci_proxy_good` | no |
| `nm_like_candidate` | no |
| `prl_like_strong_candidate` | no |

Vertex-level correlations:

| Correlation | Value |
|---|---:|
| `corr(Theta_i, abs(B_i))` | `0.5730815442` |
| `corr(Theta_i, abs(du_i))` | `0.5326892507` |
| `corr(Theta_i, abs(K_lin_i))` | `0.2884886315` |
| `corr(abs(B_i), abs(du_i))` | `0.8792203134` |

Interpretation: `Theta_NM` is not already low. The Ricci/linear proxy is not merely overconstraining a good NM configuration in this run.

## Theta-Aware Coordinate Optimization

Command:

```bash
./inc_ricci_v13 --nx 16 --ny 16 \
  --radius_seed 101 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --coord_mode soft_relax \
  --target_phi 0.84 \
  --relax_steps 30000 \
  --use_voropp 1 \
  --compute_nm_theta 1 \
  --theta_coord_opt 1 \
  --theta_opt_steps 2000 \
  --theta_opt_trials 12 \
  --theta_opt_move_scale 0.02 \
  --theta_opt_temperature 0.01 \
  --out theta_opt_v13_16_delta025_ratio5_phi084 \
  --progress_every 100
```

| Diagnostic | Initial | Final |
|---|---:|---:|
| `Theta_mean_abs` | `0.118152` | `0.117376` |
| `Theta_rms` | `0.180370` | `0.179703` |
| `Theta_max_abs` | `0.796788` | `0.796188` |
| `B_rms` | `0.528627` | `0.519804` |
| `S_pred` | `129.9249` | `129.4004` |
| `du_rms` | `0.140375` | `0.139132` |
| `du_max` | `0.615459` | `0.616352` |
| `S_theta` | `2.74197` | `2.73077` |

Theta-aware coordinate optimization slightly reduces `Theta_NM`, `B_rms`, `S_pred`, and `du_rms`, but the change is small. It does not create a divergence where theta becomes low while Ricci proxy remains high.

## N=256 Theta Multistart

Command:

```bash
./inc_ricci_v13 --nx 16 --ny 16 \
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
  --compute_nm_theta 1 \
  --theta_multistart 1 \
  --threads 16 \
  --multistart_count 32 \
  --theta_multistart_out theta_multistart_v13_16.csv \
  --progress_every 0
```

| Statistic | Value |
|---|---:|
| rows | 32 |
| `nm_theta_low` | `0 / 32` |
| `ricci_proxy_good` | `0 / 32` |
| `nm_like_candidate` | `0 / 32` |
| `prl_like_strong_candidate` | `0 / 32` |
| best `Theta_mean_abs` | `0.1105432785` |
| median `Theta_mean_abs` | `0.1199127668` |
| best `du_rms` | `0.1312777890` |
| median `du_rms` | `0.1452826070` |
| best `S_pred` | `81.95102676` |

Configuration-level correlations:

| Correlation | Value |
|---|---:|
| `corr(Theta_mean_abs, du_rms)` | `0.6093` |
| `corr(Theta_mean_abs, S_pred)` | `0.3339` |
| `corr(Theta_mean_abs, B_rms)` | `0.5406` |
| `corr(Theta_mean_abs, U_overlap)` | `0.3526` |

Interpretation: lower theta tends to accompany lower `du_rms` and lower angle-budget defect. The two diagnostics are not identical, but they are not decoupled.

## N=576 Theta Multistart

Command:

```bash
./inc_ricci_v13 --nx 24 --ny 24 \
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
  --compute_nm_theta 1 \
  --theta_multistart 1 \
  --threads 16 \
  --multistart_count 16 \
  --theta_multistart_out theta_multistart_v13_24.csv \
  --progress_every 0
```

| Statistic | Value |
|---|---:|
| rows | 16 |
| `nm_theta_low` | `0 / 16` |
| `ricci_proxy_good` | `0 / 16` |
| `nm_like_candidate` | `0 / 16` |
| `prl_like_strong_candidate` | `0 / 16` |
| best `Theta_mean_abs` | `0.1150546792` |
| median `Theta_mean_abs` | `0.1186609999` |
| best `du_rms` | `0.1384877873` |
| median `du_rms` | `0.1440729106` |
| best `S_pred` | `97.13637540` |

Configuration-level correlations:

| Correlation | Value |
|---|---:|
| `corr(Theta_mean_abs, du_rms)` | `0.6324` |
| `corr(Theta_mean_abs, S_pred)` | `0.3154` |
| `corr(Theta_mean_abs, B_rms)` | `0.3811` |
| `corr(Theta_mean_abs, U_overlap)` | `-0.0390` |

The N=576 trend matches N=256: theta and Ricci response are moderately aligned, while overlap energy is not a reliable proxy for theta.

## Theta Phi Scan, N=256

Command used `phi_list=0.76,0.78,0.80,0.82,0.84,0.86`.

| phi | `Theta_mean_abs` | `du_rms` | `du_max` | `B_rms` | `S_pred` | max overlap | theta low |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.76 | `0.1971` | `0.2014` | `0.7204` | `0.7296` | `204.6579` | `5.89e-08` | no |
| 0.78 | `0.1839` | `0.2092` | `0.6190` | `0.7852` | `184.4869` | `8.69e-08` | no |
| 0.80 | `0.1613` | `0.1795` | `0.5567` | `0.6189` | `144.0707` | `9.33e-08` | no |
| 0.82 | `0.1432` | `0.1665` | `0.5672` | `0.6191` | `136.4960` | `9.72e-08` | no |
| 0.84 | `0.1216` | `0.1444` | `0.4980` | `0.5380` | `106.4325` | `0.00387` | no |
| 0.86 | `0.1068` | `0.1423` | `0.5082` | `0.5229` | `107.0140` | `0.04796` | no |

Higher density improves both `Theta_NM` and `du_rms`, but also increases residual overlaps. None of the tested densities reaches `Theta_mean_abs < 0.05`.

## Conclusion

v13 resolves the ambiguity from v12:

```text
Theta_NM low while Ricci proxy high: no
Theta_NM high while Ricci proxy high: yes
Theta_NM correlated with du_rms/S_pred: moderately yes
```

The Ricci/linear-response proxy is not simply too strict relative to NM-style steric angle order. In the tested configurations, the true coordinate-level steric angle residual is also high.

Current interpretation:

- `Theta_NM` and the Ricci proxy are measuring related but not identical defects.
- `Theta_NM` is moderately correlated with `du_rms` across both N=256 and N=576 multistarts.
- The best `Theta_mean_abs` remains around `0.11`, above the moderate threshold `0.05`.
- Direct small-step theta-aware coordinate optimization barely improves theta.
- Overlap relaxation alone is not predictive of theta quality.

Next direction:

The next useful step should directly optimize topology/coordinates for `Theta_NM`, likely with larger topology-changing coordinate moves, particle swaps, annealing, or an explicit theta/local-motif objective. The result does not support abandoning the Ricci proxy; it supports adding true theta as a parallel optimization target.
