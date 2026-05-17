# INC-Ricci v17 defect-channel operator atlas report

## Build And Tests

Build command used under WSL:

```bash
g++ -O3 -std=c++17 inc_ricci_v17.cpp -o inc_ricci_v17 \
  -I/home/tianyu/miniconda3/envs/analyse/include \
  -I/home/tianyu/miniconda3/envs/analyse/include/voro++ \
  -L/home/tianyu/miniconda3/envs/analyse/lib \
  -lvoro++ \
  -Wl,-rpath,/home/tianyu/miniconda3/envs/analyse/lib \
  -pthread
```

Test command:

```bash
./inc_ricci_v17 --test
```

Result: all v17 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Hex baseline | PASS | `Theta_sum=2.41e-16`, `Theta_shape=3.70e-16`, `du_rms=2.76e-22` |
| Free Ricci on radius-noised hex | PASS | `Theta_sum 2.41e-02 -> 7.81e-12` |
| Shape descent small | PASS | trajectory written |
| Small operator atlas | PASS | `operator_atlas.csv`, `operator_statistics.csv` produced |
| Fine compression small | PASS | CSV produced |

## v17 Implementation Summary

`inc_ricci_v17.cpp` preserves v16 and adds:

- `operator_atlas.csv` records with before/after defect-plane displacements;
- operator families: sum-channel, shape-channel, coupled, topology;
- free, bounded, and regularized Ricci diagnostic operators;
- shape descent, fixed-topology shape descent, box aspect strain, compression scan, post-voro surgery, radius-label swap, and combined-cycle operators;
- threaded operator multistart;
- aggregate `operator_statistics.csv` / `.txt`;
- `plot_operator_vectors.py`;
- fine compression scan;
- local field correlation output.

## Single N=256 Full Operator Atlas

Command:

```bash
./inc_ricci_v17 --nx 16 --ny 16 --seed 1 --radius_seed 101 \
  --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 \
  --coord_mode soft_relax --target_phi 0.84 --relax_steps 20000 \
  --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 \
  --operator_atlas 1 --operator_set full \
  --out operator_atlas_v17_16_full --progress_every 0
```

Baseline point:

```text
Theta_sum  ~= 0.0666
Theta_shape ~= 0.1175
du_rms     ~= 0.1317
```

Selected displacement vectors:

| Operator | `Delta Theta_sum` | `Delta Theta_shape` | `Delta du_rms` | Main direction |
|---|---:|---:|---:|---|
| `free_ricci` | `-0.0666` | `-0.0025` | `-0.1317` | strong left, slight down |
| `bounded_ricci_du010` | `-0.0144` | `-0.0011` | `-0.0239` | controlled left |
| `bounded_ricci_du250` | `-0.0349` | `-0.0022` | `-0.0570` | stronger left, more overlap/radius change |
| `regularized_ricci_lsort100` | `-0.0652` | `-0.0026` | `-0.1149` | strong left |
| `shape_descent` | `0.0000` | `-0.00043` | `0.0000` | weak down |
| `fixed_topology_shape_descent` | `~0` | `-0.00037` | `0.0000` | weak pure down |
| `box_shear_aspect` | `0.0000` | `+0.000002` | `0.0000` | neutral |
| `compression_step` | `0.0000` | `~0` | `0.0000` | neutral for this seed |
| `post_voro_surgery` | `-0.00010` | `+0.0306` | `-0.00365` | left but strongly up |
| `radius_label_swap` | `-0.00004` | `-0.00001` | `-0.00014` | tiny diagonal |
| `combined_cycle` | `-0.0144` | `-0.0013` | `-0.0239` | bounded-Ricci-like diagonal |

Interpretation: hard/Ricci operators are clean sum-channel operators. Post-voro topology surgery improves hard closure but damages the NM shape channel.

## N=256 Basic Operator Multistart

Command used 96 starts and 20 threads. Result: `672` operator rows.

| Operator | mean `Delta_sum` | mean `Delta_shape` | left | down | diagonal | mean `Delta_du_rms` | overlap cost |
|---|---:|---:|---:|---:|---:|---:|---:|
| `free_ricci` | `-0.0738` | `-0.00354` | `1.00` | `1.00` | `1.00` | `-0.1468` | `+4.08` |
| `bounded_ricci_du010` | `-0.0179` | `-0.00147` | `1.00` | `1.00` | `1.00` | `-0.0294` | `+0.118` |
| `post_voro_surgery` | `-0.00900` | `+0.0295` | `0.990` | `0.00` | `0.00` | `-0.0142` | `0` |
| `shape_descent` | `-9.86e-06` | `-1.50e-05` | `0.177` | `0.594` | `0.125` | `+9.44e-06` | `+4.01e-05` |
| `compression_step` | `-3.32e-05` | `+5.64e-05` | `0.448` | `0.375` | `0.146` | `-1.03e-04` | `-0.00114` |
| `radius_label_swap` | `-8.89e-06` | `-1.99e-06` | `0.333` | `0.323` | `0.219` | `-1.34e-05` | `+2.46e-06` |
| `box_shear_aspect` | `-1.25e-06` | `-7.08e-06` | `0.073` | `0.240` | `0.052` | `+1.22e-05` | `-1.82e-06` |

No N=256 basic operator reached `low_sum_low_shape`.

## N=256 Full Operator Multistart

Command used 64 starts and 20 threads. Result: `960` operator rows.

Best directional classes:

| Operator | mean `Delta_sum` | mean `Delta_shape` | diagonal fraction | mean `Delta_du_rms` | radius-ratio factor |
|---|---:|---:|---:|---:|---:|
| `free_ricci` | `-0.0735` | `-0.00359` | `1.00` | `-0.1448` | `1.218` |
| `regularized_ricci_lsort10` | `-0.0723` | `-0.00360` | `1.00` | `-0.1322` | `1.177` |
| `regularized_ricci_lsort1000` | `-0.0684` | `-0.00361` | `1.00` | `-0.1159` | `1.118` |
| `bounded_ricci_du250` | `-0.0453` | `-0.00305` | `1.00` | `-0.0735` | `1.007` |
| `bounded_ricci_du010` | `-0.0186` | `-0.00154` | `1.00` | `-0.0305` | `0.982` |
| `combined_cycle` | `-0.0186` | `-0.00155` | `1.00` | `-0.0305` | `0.982` |
| `post_voro_surgery` | `-0.00863` | `+0.0287` | `0.00` | `-0.0133` | `1.000` |
| `shape_descent` | `+1.35e-05` | `-1.55e-05` | `0.063` | `+3.66e-07` | `1.000` |

Full-set interpretation:

- clean sum-channel: `free_ricci`, `regularized_ricci`, bounded Ricci;
- clean shape-channel: fixed-topology / coordinate shape descent, but the effect is very small;
- hard-improving but shape-damaging: `post_voro_surgery`;
- best practical diagonal among distribution-preserving-ish moves: bounded Ricci, especially `du=0.10` to `0.25`;
- no operator reached the `low_sum_low_shape` region.

## Larger N=576 Basic Operator Multistart

Command used 32 starts and 20 threads. Result: `224` rows.

| Operator | mean `Delta_sum` | mean `Delta_shape` | left | down | diagonal | mean `Delta_du_rms` |
|---|---:|---:|---:|---:|---:|---:|
| `free_ricci` | `-0.0746` | `-0.00353` | `1.00` | `1.00` | `1.00` | `-0.1475` |
| `bounded_ricci_du010` | `-0.0173` | `-0.00146` | `1.00` | `1.00` | `1.00` | `-0.0284` |
| `post_voro_surgery` | `-0.00781` | `+0.0218` | `1.00` | `0.00` | `0.00` | `-0.0106` |
| `shape_descent` | `+1.24e-05` | `-1.69e-05` | `0.125` | `0.750` | `0.063` | `-1.13e-05` |
| `compression_step` | `+0.000171` | `-1.88e-05` | `0.344` | `0.469` | `0.125` | `+0.000199` |
| `radius_label_swap` | `-1.28e-05` | `-3.11e-06` | `0.250` | `0.250` | `0.125` | `-3.86e-05` |

The operator directions persist at N=576. The expensive part is still finite-difference linear response for the after-states, but the threaded run completed successfully.

## Fine Compression Scan

Command used 12 phi values and 5 seeds.

Mean values by phi:

| phi | mean `Theta_NM` | mean `Theta_sum` | mean `Theta_shape` | mean `du_rms` | mean `U_overlap` |
|---:|---:|---:|---:|---:|---:|
| 0.74 | `0.2043` | `0.1084` | `0.1829` | `0.2171` | `1.73e-14` |
| 0.80 | `0.1607` | `0.0903` | `0.1488` | `0.1826` | `1.05e-13` |
| 0.84 | `0.1183` | `0.0731` | `0.1162` | `0.1450` | `3.48e-04` |
| 0.86 | `0.1091` | `0.0709` | `0.1089` | `0.1396` | `0.0601` |
| 0.88 | `0.1029` | `0.0638` | `0.1045` | `0.1290` | `0.3536` |
| 0.90 | `0.0963` | `0.0636` | `0.0969` | `0.1272` | `0.8254` |

Per-sample optima from `fine_compression_summary.txt`:

```text
phi minimizing Theta_NM:        0.90
phi minimizing Theta_sum:       0.88
phi minimizing Theta_shape:     0.90
phi minimizing overlap-weighted score: 0.84
```

Compression gives the clearest diagonal movement toward the origin, but the high-phi improvement is paid for by substantial overlap. With the default overlap-weighted score, `phi=0.84` remains the best compromise.

## Field Correlations

Representative field-correlation run:

```bash
./inc_ricci_v17 --nx 16 --ny 16 --seed 1 --radius_seed 101 \
  --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 \
  --coord_mode soft_relax --target_phi 0.84 --relax_steps 20000 \
  --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 \
  --compute_field_correlations 1 --corr_bins 50 \
  --out field_corr_v17_16_phi084 --progress_every 0
```

Output: `local_field_correlations_v17_16_phi084.csv`.

Summary:

```text
nonempty radial bins: 47
C_cross min: -0.00221
C_cross max:  0.000756
```

`Theta_sum` and `Theta_shape` defects have weak mixed spatial cross-correlation. The first populated neighbor shell is mildly positive, but the sign is not uniform over distance, so the two defect fields are not simply the same local field.

## Conclusion

v17 identifies operator directions in the defect plane:

```text
Theta_sum   = radius-topology mismatch
Theta_shape = coordinate angular disorder
```

Operator taxonomy:

- clean sum-channel operators: `free_ricci`, `regularized_ricci`, `bounded_ricci`;
- weak shape-channel operators: `shape_descent`, `fixed_topology_shape_descent`;
- coupled diagonal operator: compression, but only with overlap cost;
- hard-improving / soft-damaging operator: `post_voro_surgery`;
- diagnostic but weak operator: `radius_label_swap`;
- no tested operator reaches `low_sum_low_shape`.

Interpretation:

The Ricci hard proxy is still best understood as a `Theta_sum` detector. It can be improved reliably by radius response or topology surgery, but topology surgery can worsen `Theta_shape`. The missing geometric mechanism for true INC-like states is an operation that lowers coordinate angular disorder at comparable strength to how Ricci lowers radius-topology mismatch, without introducing overlap or radius escape.
