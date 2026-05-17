# INC-Ricci v11 voro++ coordinate-topology report

## Build And voro++ Detection

Detected WSL conda environment components:

| Item | Path |
|---|---|
| `voro++` executable | `/home/tianyu/miniconda3/envs/analyse/bin/voro++` |
| headers | `/home/tianyu/miniconda3/envs/analyse/include/voro++/voro++.hh` |
| library | `/home/tianyu/miniconda3/envs/analyse/lib/libvoro++.a` |

Build command used:

```bash
g++ -O3 -std=c++17 inc_ricci_v11.cpp -o inc_ricci_v11 \
  -I/home/tianyu/miniconda3/envs/analyse/include \
  -I/home/tianyu/miniconda3/envs/analyse/include/voro++ \
  -L/home/tianyu/miniconda3/envs/analyse/lib \
  -lvoro++ \
  -Wl,-rpath,/home/tianyu/miniconda3/envs/analyse/lib \
  -pthread
```

`inc_ricci_v11.cpp` links against the voro++ C++ library directly. No KNN fake result is used for the voro path.

## Tests

Command:

```bash
./inc_ricci_v11 --test
```

Result: all v11 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Triangular lattice equal radii, voro++ topology | PASS | valid torus, `E=192`, `F=128`, `B_rms=1.776357e-15`, `du_rms=2.761069e-22` |
| Poisson-like lognormal coordinates | PASS | voro++ ran, valid controlled status |
| Soft relaxation | PASS | `U_initial=4.598484`, `U_final=0.012405` |
| Radical reconstruction summary | PASS | `candidate_cliques=128`, `accepted_faces=128`, valid torus |
| v11 strategy compare | PASS | `v11_strategy_compare.csv` produced |

## Implementation Summary

v11 preserves the v9/v10 linear-response diagnostics and replaces the failed KNN coordinate topology path with a voro++ radical Voronoi workflow:

- 3D thin-slab `container_poly` with periodic x/y and non-periodic z;
- radical Voronoi neighbor extraction from `voronoicell_neighbor::neighbors()`;
- periodic canonical edge shifts;
- dual triangle reconstruction from 3-cliques;
- 2D radical empty-circle validation;
- torus triangulation validation with `E=3N`, `F=2N`, Euler characteristic 0, and every edge incident to two faces;
- linear-response predictor on valid voro radical triangulations;
- optional post-voro local linear surgery, reported as no longer exact radical Delaunay.

## N=256 Single voro++ Run

Command:

```bash
./inc_ricci_v11 --nx 16 --ny 16 \
  --radius_seed 101 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --coord_mode soft_relax \
  --target_phi 0.84 \
  --relax_steps 50000 \
  --use_voropp 1 \
  --topology_from_coords voro_radical \
  --lin_lambda 1e-2 \
  --out voro_v11_16_delta025_ratio5_phi084 \
  --progress_every 0
```

Coordinate relaxation:

| Diagnostic | Value |
|---|---:|
| `U_overlap_initial` | `18.1089029494` |
| `U_overlap_final` | `5.3316467103e-05` |
| `max_overlap_final` | `1.3512276743e-03` |

voro++ graph and dual triangulation:

| Diagnostic | Value |
|---|---:|
| `N,E,F` | `256,768,512` |
| `voro_neighbor_edges` | `768` |
| `candidate_cliques` | `512` |
| `accepted_triangles` | `512` |
| Euler characteristic | `0` |
| edge incidence bad count | `0` |
| valid triangulated torus | yes |
| degree range | `4..9` |
| disconnected components | `1` |

Linear-response diagnostics:

| Diagnostic | Value |
|---|---:|
| `B_rms` | `0.5137172207` |
| `B_max` | `1.8962043324` |
| `S_pred` | `107.8701770` |
| `rms_K_lin` | `8.2876802769e-04` |
| `max_K_lin` | `3.0657922189e-03` |
| `du_rms` | `0.1388875210` |
| `du_max` | `0.5226173840` |
| degree-radius corr | `0.7954230058` |
| `linear_predicted_good` | no |

Interpretation: the voro++ radical triangulation is valid and geometrically credible. It lowers the required response compared with random abstract topology, but `du_rms` and especially `du_max` remain above the strict physical thresholds.

## Strategy Comparison, N=256

Command:

```bash
./inc_ricci_v11 --nx 16 --ny 16 \
  --seed 1 \
  --radius_seed 101 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --coord_mode soft_relax \
  --target_phi 0.84 \
  --relax_steps 50000 \
  --use_voropp 1 \
  --v11_strategy_compare 1 \
  --out strategy_v11_16_delta025_ratio5 \
  --progress_every 0
```

| Strategy | Status | `S_pred` | `du_rms` | `du_max` | `B_rms` | degree corr | predicted good |
|---|---|---:|---:|---:|---:|---:|---:|
| `abstract_v9_best` | ok | `122.6557` | `0.1756` | `0.4654` | `0.4879` | `0.8082` | no |
| `voro_poisson` | ok | `301.2367` | `0.2608` | `0.8342` | `0.9149` | `0.5974` | no |
| `voro_soft` | ok | `107.8702` | `0.1389` | `0.5226` | `0.5137` | `0.7954` | no |
| `voro_soft_post_surgery` | not radical after surgery | `59.5487` | `0.1231` | `0.2782` | `0.4595` | `0.8392` | no |
| `triangular_lattice_equal` | ok | `8.486e-13` | `2.762e-22` | `2.852e-22` | `1.776e-15` | `1.0000` | yes |

voro++ soft-relaxed radical topology beats `abstract_v9_best` in `S_pred` and `du_rms`, but not in `du_max`. Light post-voro surgery gives the best score and response, but the topology is no longer an exact radical Delaunay result.

## N=256 voro++ Multistart

Command:

```bash
./inc_ricci_v11 --nx 16 --ny 16 \
  --seed 1000 \
  --radius_seed 101 \
  --same_radii_all_starts 1 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --coord_mode soft_relax \
  --target_phi 0.84 \
  --relax_steps 30000 \
  --use_voropp 1 \
  --voro_multistart 1 \
  --threads 16 \
  --multistart_count 32 \
  --voro_multistart_out voro_multistart_v11_16_delta025_ratio5.csv \
  --progress_every 0
```

| Statistic | Value |
|---|---:|
| rows | 32 |
| valid torus rate | `32 / 32` |
| `linear_predicted_good` | `0 / 32` |
| best `du_rms` | `0.1298103177` |
| median `du_rms` | `0.1452121718` |
| best `du_max` | `0.5146946675` |
| median `du_max` | `0.5113018421` |
| fraction beating v9 best `du_rms ~= 0.133565` | `2 / 32` |
| best `S_pred` | `80.55072432` |

The voro++ reconstruction is robust in this N=256 scan. It sometimes beats the v9 best in `du_rms`, but none of the starts reach the strict predictor thresholds.

## Lower Density Check, N=256, Phi 0.78

Command used `target_phi=0.78`, 16 starts.

| Statistic | Value |
|---|---:|
| valid torus rate | `16 / 16` |
| `linear_predicted_good` | `0 / 16` |
| best `du_rms` | `0.1856178442` |
| median `du_rms` | `0.1951375248` |
| best `du_max` | `0.7427061984` |
| median `du_max` | `0.6699430756` |

Lower density was worse for the linear-response target in this implementation.

## Larger N=576 Multistart

Command:

```bash
./inc_ricci_v11 --nx 24 --ny 24 \
  --seed 2000 \
  --radius_seed 202 \
  --same_radii_all_starts 1 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --coord_mode soft_relax \
  --target_phi 0.84 \
  --relax_steps 50000 \
  --use_voropp 1 \
  --voro_multistart 1 \
  --threads 16 \
  --multistart_count 16 \
  --voro_multistart_out voro_multistart_v11_24_delta025_ratio5.csv \
  --progress_every 0
```

| Statistic | Value |
|---|---:|
| rows | 16 |
| valid torus rate | `16 / 16` |
| `linear_predicted_good` | `0 / 16` |
| best `du_rms` | `0.1372715439` |
| median `du_rms` | `0.1457434808` |
| best `du_max` | `0.5343578813` |
| median `du_max` | `0.5400866466` |
| runtime | about `172.5 s` |

The N=576 trend is consistent with N=256: topology reconstruction remains valid, but the required perturbation is still too large.

## Conclusion

v11 fixes the primary v10 bottleneck. The voro++ radical Voronoi neighbor graph plus weighted empty-circle dual reconstruction now produces valid periodic triangulated torus topologies for the tested N=256 and N=576 systems.

The result is scientifically useful but not a physical success:

```text
voro++ radical topology valid: yes
beats abstract random topology: yes
beats abstract v9 best in du_rms: sometimes
reaches strict small-du physical threshold: no
post-voro surgery improves response: yes, but topology is no longer exact radical Delaunay
```

Interpretation:

- Real-space radical topology is a real missing ingredient compared with v10 KNN/fallback topology.
- It improves `S_pred` and `du_rms`, especially after soft coordinate relaxation.
- It does not by itself solve the physical INC-like problem at `target_delta=0.25`, `target_radius_ratio=5`.
- The remaining obstruction is the required radius perturbation, especially `du_max`.
- The next useful direction is either better coordinate relaxation toward true jammed/radical structures or constrained post-voro topology optimization that preserves geometric validity.
