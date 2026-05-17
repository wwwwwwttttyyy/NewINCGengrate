# INC-Ricci v15 defect-decomposition report

## Build And Tests

Build command used under WSL:

```bash
g++ -O3 -std=c++17 inc_ricci_v15.cpp -o inc_ricci_v15 \
  -I/home/tianyu/miniconda3/envs/analyse/include \
  -I/home/tianyu/miniconda3/envs/analyse/include/voro++ \
  -L/home/tianyu/miniconda3/envs/analyse/lib \
  -lvoro++ \
  -Wl,-rpath,/home/tianyu/miniconda3/envs/analyse/lib \
  -pthread
```

Test command:

```bash
./inc_ricci_v15 --test
```

Result: all v15 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Hex decomposition baseline | PASS | `Theta=3.47e-16`, `sum=2.41e-16`, `shape=3.70e-16` |
| Coordinate-noise perturbation | PASS | `sum=2.61e-16`, `shape=2.90e-02` |
| Radius-noise perturbation | PASS | `sum=2.13e-02`, `B_rms=0.1679`, `du_rms=0.0450` |
| Soft-relaxed local defects/regression | PASS | `R2=0.786`, dominant predictor `actual_angle_std` |
| Topology contribution compare | PASS | 4 rows written |

## v15 Implementation Summary

`inc_ricci_v15.cpp` preserves v14 and adds:

- `local_defects.csv` with per-particle steric, topology, radius-environment, overlap/gap, and angular-environment diagnostics;
- local Pearson correlations and ridge OLS regression;
- theta decomposition into:
  - `Theta_sum_component`: total angle-budget/radius-topology mismatch;
  - `Theta_shape_component`: coordinate angular redistribution disorder;
- perturbation decomposition CSV;
- density decomposition scan;
- topology contribution comparison;
- radius-escape decomposition summary;
- automatic labels: `coordinate_dominated`, `radius_topology_dominated`, or `mixed`.

## Single N=256 Decomposition

Command:

```bash
./inc_ricci_v15 --nx 16 --ny 16 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 30000 --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 --out decomp_v15_16_delta025_phi084 --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| `Theta_NM` | `0.1181521953` |
| `Theta_sum_component` | `0.0718377065` |
| `Theta_shape_component` | `0.1167310825` |
| dominant label | `mixed` |
| `B_rms` | `0.5286270607` |
| `du_rms` | `0.1403753400` |
| `du_max` | `0.6154590162` |

Local regression:

| Field | Value |
|---|---:|
| `R2` | `0.7690946278` |
| dominant predictor | `actual_angle_std` |
| hard proxy validity | `moderate` |

Selected local correlations:

| Correlation | Value |
|---|---:|
| `corr(Theta_i, abs(B_i))` | `0.5731` |
| `corr(Theta_i, abs(du_i))` | `0.5327` |
| `corr(Theta_i, actual_angle_std_i)` | `0.7232` |
| `corr(Theta_i, ideal_angle_std_i)` | `-0.0480` |
| `corr(Theta_i, local_free_gap_min_i)` | `0.1260` |

Interpretation: ordinary soft-relaxed voro topology is not purely coordinate-dominated or purely radius-topology-dominated. The shape component is larger, but the sum component remains substantial.

## Perturbation Decomposition

Coordinate noise from hex baseline:

| amp | `Theta_NM` | sum | shape | `B_rms` | `du_rms` | regime |
|---:|---:|---:|---:|---:|---:|---|
| 0 | `3.47e-16` | `2.41e-16` | `3.70e-16` | `1.78e-15` | `2.76e-22` | `soft_low_hard_low` |
| 0.02 | `0.01367` | `2.59e-16` | `0.01367` | `1.78e-15` | `2.76e-22` | `soft_low_hard_low` |
| 0.05 | `0.03565` | `2.58e-16` | `0.03565` | `1.78e-15` | `2.76e-22` | `soft_low_hard_low` |
| 0.10 | `0.07015` | `2.53e-16` | `0.07015` | `1.78e-15` | `2.76e-22` | `soft_high_hard_low` |
| 0.20 | `0.13889` | `2.58e-16` | `0.13889` | `1.78e-15` | `2.76e-22` | `soft_high_hard_low` |

Coordinate noise is almost a pure shape-component perturbation. It can make the sample soft-bad while the hard Ricci proxy stays perfect.

Radius noise from hex baseline:

| amp | `Theta_NM` | sum | shape | `B_rms` | `du_rms` | `du_max` | regime |
|---:|---:|---:|---:|---:|---:|---:|---|
| 0.02 | `0.01191` | `0.01057` | `0.00536` | `0.07844` | `0.02013` | `0.06335` | `soft_low_hard_low` |
| 0.05 | `0.03069` | `0.02744` | `0.01448` | `0.20281` | `0.05376` | `0.18722` | `soft_low_hard_high` |
| 0.10 | `0.06194` | `0.05338` | `0.02893` | `0.41823` | `0.10804` | `0.36448` | `soft_high_hard_high` |
| 0.20 | `0.10295` | `0.08888` | `0.05088` | `0.68182` | `0.19091` | `0.64276` | `soft_high_hard_high` |

Radius noise primarily increases the sum component and the hard metrics. This is the clean radius-topology mismatch channel.

## Density Decomposition Scan

| phi | `Theta_NM` | sum | shape | `B_rms` | `du_rms` | `du_max` | max overlap |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.76 | `0.1971` | `0.0975` | `0.1812` | `0.7296` | `0.2014` | `0.7204` | `5.89e-08` |
| 0.80 | `0.1613` | `0.0841` | `0.1557` | `0.6189` | `0.1795` | `0.5567` | `9.33e-08` |
| 0.84 | `0.1216` | `0.0744` | `0.1203` | `0.5380` | `0.1444` | `0.4980` | `3.87e-03` |
| 0.86 | `0.1068` | `0.0731` | `0.1068` | `0.5229` | `0.1423` | `0.5082` | `4.80e-02` |
| 0.88 | `0.0987` | `0.0612` | `0.1017` | `0.4509` | `0.1207` | `0.3311` | `9.18e-02` |

Increasing density reduces both components, with a stronger visible reduction in the shape component. The tradeoff is increasing residual overlap.

## Topology Contribution

| Topology | `Theta_ring` | sum | shape | `B_rms` | `S_pred` | `du_rms` | `du_max` |
|---|---:|---:|---:|---:|---:|---:|---:|
| radical voro | `0.1182` | `0.0718` | `0.1167` | `0.5286` | `129.925` | `0.1404` | `0.6155` |
| post-surgery nonradical | `0.1442` | `0.0654` | `0.1507` | `0.4719` | `63.287` | `0.1270` | `0.2923` |
| random flip | `0.6797` | `0.2393` | `0.6337` | `1.7218` | `2181.21` | `0.5090` | `2.8229` |
| abstract degree-angle | `0.6690` | `0.1150` | `0.6533` | `0.6503` | `138.414` | `0.1834` | `0.5153` |

The radical topology strongly improves both soft and hard metrics relative to random or abstract topology at fixed coordinates. Light post-surgery improves hard closure and the sum component, but worsens shape/theta and is no longer radical Delaunay.

## Radius Escape Decomposition

| Diagnostic | Initial | Final Free Ricci |
|---|---:|---:|
| `Theta_NM` | `0.117236` | `0.112137` |
| `Theta_sum_component` | `0.066591` | `6.83e-12` |
| `Theta_shape_component` | `0.117539` | `0.115018` |
| `B_rms` | `0.488451` | `4.60e-11` |
| `du_rms` | `0.131671` | `2.15e-10` |
| radius ratio | `3.96077` | `5.33552` |

Free Ricci closure almost entirely removes the sum component while leaving the shape component nearly unchanged. Radius relaxation is therefore directly relaxing radius-topology angle-budget mismatch, not coordinate angular disorder.

## Conclusion

v15 decomposes the defect into useful channels:

```text
coordinate noise channel: shape component
radius noise channel: sum component + hard proxy
ordinary soft-relaxed voro samples: mixed, shape larger but sum still substantial
free Ricci escape: removes sum component, not shape component
```

Interpretation:

- The missing ingredient is not only radius-topology matching and not only coordinate angular ordering.
- For the tested `Delta=0.25`, ratio `5` soft-relaxed radical samples, both components matter.
- Hard Ricci closure mostly targets the sum/radius-topology channel.
- NM-style `Theta_NM` also sees coordinate shape disorder, which can be invisible to the hard proxy if topology and radii remain ideal.
- A future constructive search should optimize both: radius/topology compatibility and coordinate angular regularity.
