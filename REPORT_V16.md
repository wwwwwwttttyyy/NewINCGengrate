# INC-Ricci v16 defect-plane atlas report

## Build And Tests

Build command used under WSL:

```bash
g++ -O3 -std=c++17 inc_ricci_v16.cpp -o inc_ricci_v16 \
  -I/home/tianyu/miniconda3/envs/analyse/include \
  -I/home/tianyu/miniconda3/envs/analyse/include/voro++ \
  -L/home/tianyu/miniconda3/envs/analyse/lib \
  -lvoro++ \
  -Wl,-rpath,/home/tianyu/miniconda3/envs/analyse/lib \
  -pthread
```

Test command:

```bash
./inc_ricci_v16 --test
```

Result: all v16 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Hex baseline | PASS | `Theta_sum=2.41e-16`, `Theta_shape=3.70e-16`, `du_rms=2.76e-22` |
| Hex coord noise | PASS | `Theta_shape=2.90e-02 > Theta_sum=2.61e-16` |
| Hex radius noise | PASS | `Theta_sum=2.13e-02`, `B_rms=1.68e-01`, `du_rms=4.50e-02` |
| Small atlas | PASS | `defect_plane_atlas_v16_test.csv` produced |
| Plot script | PASS | `plot_defect_plane.py` generated |

## v16 Implementation Summary

`inc_ricci_v16.cpp` preserves the v15 decomposition pipeline and adds:

- unified `DefectPlaneRecord`;
- `defect_plane_atlas.csv` style export;
- defect-plane regime classification from `Theta_sum` and `Theta_shape`;
- soft/hard regime classification from `Theta_NM`, `B_rms`, `du_rms`, and `du_max`;
- atlas classes for hex baseline, coordinate noise, radius noise, combined noise, phi scan, radius escape endpoints, and controlled placeholder rows for unavailable nonradical classes;
- trajectory CSVs for coordinate perturbation, radius perturbation, combined perturbation, phi scan, and radius escape;
- atlas-level correlations and two-predictor regressions;
- `plot_defect_plane.py`.

## N=256 Defect-Plane Atlas

Command:

```bash
./inc_ricci_v16 --nx 16 --ny 16 \
  --radius_seed 101 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --build_defect_atlas 1 \
  --coord_mode soft_relax \
  --relax_steps 20000 \
  --use_voropp 1 \
  --compute_nm_theta 1 \
  --decompose_theta 1 \
  --defect_plane_out defect_plane_atlas_v16_16.csv \
  --out defect_plane_v16_16 \
  --progress_every 0
```

Result: `33` atlas records.

| Defect-plane regime | Count |
|---|---:|
| `low_sum_low_shape` | 17 |
| `low_sum_high_shape` | 2 |
| `high_sum_low_shape` | 0 |
| `high_sum_high_shape` | 12 |
| `unknown` | 2 |

| Soft-hard regime | Count |
|---|---:|
| `soft_low_hard_low` | 12 |
| `soft_low_hard_high` | 4 |
| `soft_high_hard_low` | 2 |
| `soft_high_hard_high` | 13 |
| `unknown` | 2 |

Sample-class locations:

| Class | Region summary | Notes |
|---|---|---|
| `hex_baseline` | `low_sum_low_shape` | zero baseline fixed |
| `hex_coord_noise` | `low_sum_low_shape` to `low_sum_high_shape` | vertical motion in defect plane |
| `hex_radius_noise` | mostly `low_sum_low_shape`, high amplitude `high_sum_high_shape` | hard metrics rise with `Theta_sum` |
| `hex_combined_noise` | mixed, reaches `high_sum_high_shape` | both channels activated |
| `voro_soft_phi_scan` | all `high_sum_high_shape` | ordinary soft-relaxed radical samples remain far from either ideal |
| `post_surgery_nonradical`, `abstract_v9_like` | `unknown` | intentionally labeled nonradical / placeholder |

## Controlled Perturbations

Hex coordinate noise moves almost purely in the shape direction:

| amp | `Theta_NM` | `Theta_sum` | `Theta_shape` | `B_rms` | `du_rms` |
|---:|---:|---:|---:|---:|---:|
| 0 | `3.47e-16` | `2.41e-16` | `3.70e-16` | `1.78e-15` | `2.76e-22` |
| 0.05 | `0.0345` | `2.49e-16` | `0.0345` | `1.78e-15` | `2.76e-22` |
| 0.10 | `0.0679` | `2.57e-16` | `0.0679` | `1.78e-15` | `2.76e-22` |
| 0.20 | `0.1370` | `2.55e-16` | `0.1370` | `1.78e-15` | `2.76e-22` |

Hex radius noise moves mainly in the sum / hard direction:

| amp | `Theta_NM` | `Theta_sum` | `Theta_shape` | `B_rms` | `du_rms` | `du_max` |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | `3.47e-16` | `2.41e-16` | `3.70e-16` | `1.78e-15` | `2.76e-22` | `2.85e-22` |
| 0.05 | `0.0268` | `0.0235` | `0.0133` | `0.1728` | `0.0478` | `0.1379` |
| 0.10 | `0.0529` | `0.0457` | `0.0259` | `0.3442` | `0.0954` | `0.2715` |
| 0.20 | `0.1045` | `0.0897` | `0.0532` | `0.6979` | `0.1999` | `0.6643` |

This confirms the intended split: coordinate angular disorder is visible to `Theta_NM` but largely invisible to the Ricci hard proxy when topology and radii remain compatible.

## Phi Scan

For N=256, increasing density from `phi=0.76` to `0.86` moves the voro-soft samples toward the origin in both channels, then `phi=0.88` worsens the sum channel while overlaps increase.

| phi | `Theta_NM` | `Theta_sum` | `Theta_shape` | `du_rms` | `du_max` | `max_overlap` |
|---:|---:|---:|---:|---:|---:|---:|
| 0.76 | `0.1919` | `0.0964` | `0.1712` | `0.1931` | `0.6507` | `5.76e-08` |
| 0.82 | `0.1419` | `0.0798` | `0.1380` | `0.1597` | `0.5180` | `1.23e-07` |
| 0.84 | `0.1172` | `0.0666` | `0.1175` | `0.1317` | `0.5227` | `0.00350` |
| 0.86 | `0.1049` | `0.0642` | `0.1060` | `0.1284` | `0.3432` | `0.0481` |
| 0.88 | `0.1054` | `0.0755` | `0.1070` | `0.1434` | `0.4906` | `0.0903` |

Improvement is diagonal through `Theta_sum` and `Theta_shape`, but high density introduces overlap.

## Radius Escape

For the N=256 voro-soft topology:

| Stage | `B_rms` | `du_rms` | radius ratio | note |
|---|---:|---:|---:|---|
| initial | `0.4885` | `0.1317` | `3.9608` | `high_sum_high_shape` |
| free Ricci endpoint | `4.60e-11` | n/a | `5.3355` | curvature closed without catastrophic radius escape |

The atlas records the escape hard trajectory endpoints. The fixed coordinate-level `Theta_shape` remains high; free Ricci removes the curvature/hard mismatch but does not repair coordinate angular disorder.

## Correlations

N=256 atlas correlations:

| Quantity | Value |
|---|---:|
| `corr(Theta_NM, Theta_sum)` | `0.8875` |
| `corr(Theta_NM, Theta_shape)` | `0.9844` |
| `corr(Theta_NM, du_rms)` | `0.8802` |
| `corr(Theta_sum, B_rms)` | `0.9502` |
| `corr(Theta_sum, du_rms)` | `0.9983` |
| `corr(Theta_shape, du_rms)` | `0.8189` |

Regression:

```text
Theta_NM ~= 0.00248 + 0.36523*Theta_sum + 0.86075*Theta_shape
du_rms   ~= 0.00026 + 2.10348*Theta_sum - 0.04656*Theta_shape
```

The regression is the cleanest separation: `Theta_NM` depends on both channels, dominated by shape in this atlas, while `du_rms` is essentially a `Theta_sum` detector.

## N=576 Check

Command:

```bash
./inc_ricci_v16 --nx 24 --ny 24 \
  --radius_seed 202 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --build_defect_atlas 1 \
  --coord_mode soft_relax \
  --relax_steps 30000 \
  --use_voropp 1 \
  --compute_nm_theta 1 \
  --decompose_theta 1 \
  --defect_plane_out defect_plane_atlas_v16_24.csv \
  --out defect_plane_v16_24 \
  --progress_every 0
```

Result: `33` records. Regime counts match N=256: `17 low_sum_low_shape`, `2 low_sum_high_shape`, `12 high_sum_high_shape`, `2 unknown`.

Key N=576 correlations:

| Quantity | Value |
|---|---:|
| `corr(Theta_NM, Theta_sum)` | `0.8959` |
| `corr(Theta_NM, Theta_shape)` | `0.9813` |
| `corr(Theta_NM, du_rms)` | `0.9046` |
| `corr(Theta_sum, B_rms)` | `0.9396` |
| `corr(Theta_sum, du_rms)` | `0.9985` |
| `corr(Theta_shape, du_rms)` | `0.8419` |

Regression:

```text
Theta_NM ~= 0.00309 + 0.39398*Theta_sum + 0.82894*Theta_shape
du_rms   ~= 0.00054 + 1.92960*Theta_sum + 0.04393*Theta_shape
```

The larger system preserves the same qualitative separation.

## Conclusion

v16 turns the v15 decomposition into an explicit atlas:

```text
Theta_sum   = radius-topology mismatch channel
Theta_shape = coordinate angular disorder channel
```

Main result:

- coordinate noise moves vertically: `Theta_shape` rises while hard closure remains near zero;
- radius noise moves horizontally/diagonally: `Theta_sum`, `B_rms`, and `du_rms` rise together;
- ordinary voro-soft polydisperse samples occupy `high_sum_high_shape`;
- increasing density moves voro-soft samples toward the origin until overlap becomes significant;
- free Ricci closure removes the hard curvature mismatch but does not repair coordinate shape disorder.

Interpretation:

`Theta_NM` and PRL/Ricci hard closure are aligned only when the dominant defect is radius-topology mismatch. The Ricci proxy is mainly a detector of `Theta_sum`. True NM-like steric order requires both `Theta_sum` and `Theta_shape` to be small, so the remaining missing channel is coordinate angular ordering, not just hard radius-topology closure.
