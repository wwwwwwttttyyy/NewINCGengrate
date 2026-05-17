# INC-Ricci v10 coordinate-generated topology report

## Build And Tests

Build:

```bash
g++ -O3 -std=c++17 inc_ricci_v10.cpp -o inc_ricci_v10
```

Test command:

```bash
./inc_ricci_v10 --test
```

Result: all v10 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Regular `8x8`, equal radii | PASS | `B_rms=1.776357e-15`, `S_pred=8.486082e-13` |
| Poisson-like coordinate KNN | PASS | raw KNN invalid, controlled fallback valid |
| Soft relaxation | PASS | overlap energy decreased `4.598484 -> 0.012405` |
| Brute-force weighted Delaunay | PASS | controlled failure, no crash |
| Coordinate strategy compare | PASS | `coord_strategy_compare_v10_test.csv` produced |

## Implementation Summary

`inc_ricci_v10.cpp` preserves the v9 linear-response and topology tools, and adds:

- coordinate generation modes `random`, `poisson_like`, `soft_relax`, `none`;
- periodic soft-overlap relaxation;
- approximate KNN topology builder;
- brute-force weighted-Delaunay-like builder for small systems;
- coordinate-biased abstract fallback via edge flips;
- coordinate topology quality outputs;
- coordinate strategy comparison;
- coordinate multistart scan.

Important limitation: `coord_knn` is only an approximate geometric topology, not exact radical Delaunay. The brute-force path is explicitly marked as controlled failure when it does not produce a valid torus triangulation.

## Single Coordinate KNN Run

Command:

```bash
./inc_ricci_v10 --nx 8 --ny 8 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 20000 --topology_from_coords knn --knn_k 8 --alpha_radius 0.5 --out coord_knn_v10_8_delta025_ratio5 --progress_every 0
```

Coordinate relaxation:

| Diagnostic | Value |
|---|---:|
| `U_overlap_initial` | `3.992969` |
| `U_overlap_final` | `3.092776e-03` |
| `max_overlap_initial` | `0.693046` |
| `max_overlap_final` | `0.017916` |
| final overlap pair count | 132 |

Topology:

| Diagnostic | Value |
|---|---:|
| raw builder | `knn` |
| raw KNN status | invalid: edge has more than two incident faces |
| fallback used | yes, `abstract_seeded` |
| final `N,E,F` | `64,192,128` |
| Euler characteristic | `0` |
| edge/face count error | `0 / 0` |

Linear-response diagnostics after fallback:

| Diagnostic | Value |
|---|---:|
| `B_rms` | `1.419268` |
| `B_max` | `3.626214` |
| `S_pred` | `3760.4564` |
| `rms_K_lin` | `1.077405e-02` |
| `du_rms` | `0.920082` |
| `du_max` | `3.055659` |
| degree-radius corr | `0.353154` |
| linear predicted good | no |

Interpretation: coordinate relaxation successfully reduces overlaps, but the approximate KNN builder does not produce a valid torus triangulation. The fallback topology is legal but much worse than the v9 angle/degree search.

## Coordinate Strategy Compare

Command:

```bash
./inc_ricci_v10 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --compare_coord_topologies 1 --coord_mode soft_relax --target_phi 0.84 --relax_steps 20000 --topology_from_coords knn --knn_k 8 --alpha_radius 0.5 --out coord_strategy_compare_v10_8 --progress_every 0
```

| Strategy | Status | Valid torus | `S_pred` | `du_rms` | `du_max` | `B_rms` | degree corr |
|---|---|---:|---:|---:|---:|---:|---:|
| abstract_random | ok | 1 | `2422.3295` | `0.852836` | `1.993017` | `2.624738` | `-0.030071` |
| abstract_degree_angle | ok | 1 | `127.9419` | `0.187467` | `0.442990` | `0.541338` | `0.807544` |
| coord_knn | invalid edge incidence | 0 | n/a | n/a | n/a | n/a | n/a |
| coord_soft_knn | invalid edge incidence | 0 | n/a | n/a | n/a | n/a | n/a |
| coord_bruteforce | controlled failure | 0 | n/a | n/a | n/a | n/a | n/a |
| coord_seeded_flip | fallback valid | 1 | `3449.0038` | `0.951011` | `2.671243` | `1.347168` | `0.604973` |

Best topology source: `abstract_degree_angle`.

Coordinate-generated KNN did not beat the abstract v9 pipeline because it did not generate a valid triangulated torus. The coordinate-biased fallback was valid but worse.

## Brute-Force Weighted Delaunay Test

Command:

```bash
./inc_ricci_v10 --nx 8 --ny 8 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 20000 --topology_from_coords delaunay_bruteforce --delaunay_bruteforce 1 --out coord_delaunay_v10_8_delta025_ratio5 --progress_every 0
```

Result:

| Diagnostic | Value |
|---|---:|
| raw builder | `delaunay_bruteforce` |
| raw status | controlled failure: weighted complex is not a valid torus triangulation |
| fallback used | yes |
| final `N,E,F` | `64,192,128` |
| `B_rms` | `1.419268` |
| `S_pred` | `3760.4564` |
| `du_rms` | `0.920082` |
| `du_max` | `3.055659` |

Interpretation: the brute-force weighted complex did not produce a valid triangulated torus in this implementation. Since the result falls back to the same abstract-seeded topology, it is not better than KNN.

## Coordinate Multistart KNN, 8x8

Command:

```bash
./inc_ricci_v10 --nx 8 --ny 8 --seed 1000 --radius_seed 101 --same_radii_all_starts 1 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 10000 --topology_from_coords knn --knn_k 8 --alpha_radius 0.5 --coord_multistart 1 --threads 16 --multistart_count 100 --coord_multistart_out coord_multistart_knn_v10_8.csv --progress_every 0
```

| Statistic | Value |
|---|---:|
| rows | 100 |
| valid fallback torus | 100 / 100 |
| `linear_predicted_good` | 0 / 100 |
| best run | 54 |
| best `S_pred` | `923.3719` |
| best `du_rms` | `0.493635` |
| best `du_max` | `1.361116` |
| min `du_rms` | `0.493635` |
| median `du_rms` | `0.904825` |
| min `du_max` | `1.361116` |
| median `du_max` | `2.988518` |
| fraction beating v9 best `du_rms=0.133565` | 0 / 100 |
| corr(`U_overlap_final`, `S_pred`) | `0.0393` |

The coordinate soft-overlap energy is not predictive of the final linear-response score in this fallback path.

## Larger 12x12 Coordinate Multistart

Command:

```bash
./inc_ricci_v10 --nx 12 --ny 12 --seed 2000 --radius_seed 202 --same_radii_all_starts 1 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 15000 --topology_from_coords knn --knn_k 8 --alpha_radius 0.5 --coord_multistart 1 --threads 16 --multistart_count 32 --coord_multistart_out coord_multistart_knn_v10_12.csv --progress_every 0
```

| Statistic | Value |
|---|---:|
| rows | 32 |
| valid fallback torus | 32 / 32 |
| `linear_predicted_good` | 0 / 32 |
| best run | 12 |
| best `S_pred` | `2891.3025` |
| best `du_rms` | `0.842548` |
| best `du_max` | `2.551370` |
| min `du_rms` | `0.842548` |
| median `du_rms` | `1.051073` |
| min `du_max` | `2.551370` |
| median `du_max` | `3.490246` |
| corr(`U_overlap_final`, `S_pred`) | `0.1079` |

The larger system shows the same negative trend.

## Conclusion

v10 answers the coordinate-topology hypothesis in a qualified way:

- Real-space coordinate relaxation works as a coordinate generator.
- The approximate KNN builder does not produce a valid triangulated torus.
- The simple brute-force weighted-Delaunay-like builder also fails to produce a valid torus complex.
- The coordinate-biased abstract fallback produces legal triangulations, but they are worse than the v9 abstract degree+angle pipeline.

Current interpretation:

```text
coordinate-generated exact topology tested: not yet, because KNN/bruteforce builders failed validity
coordinate fallback beats v9: no
KNN validity problem: yes, primary bottleneck
missing ingredient: proper periodic radical Delaunay / regular triangulation implementation
```

The next useful step is not more Ricci tuning. It is a robust periodic radical-Delaunay/regular-triangulation builder, likely requiring either a carefully implemented computational-geometry module or an external geometry library in a later non-single-file version.

