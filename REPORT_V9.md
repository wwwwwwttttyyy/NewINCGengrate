# INC-Ricci v9 angle-budget topology-design report

## Build And Tests

Build:

```bash
g++ -O3 -std=c++17 inc_ricci_v9.cpp -o inc_ricci_v9
```

Test command:

```bash
./inc_ricci_v9 --test
```

Result: all v9 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Regular `8x8`, equal radii | PASS | `B_rms=1.776357e-15`, `S_pred=8.486082e-13` |
| Random `8x8`, fixed target radii | PASS | `B_rms=2.091344`, `B_max=7.165352` |
| Small angle precondition | PASS | `E_B` decreased `127.8086 -> 10.8188` |
| v9 strategy compare | PASS | 7 strategy rows written |
| Parallel multistart small | PASS | 4 rows written, no race issues |

## Implementation Summary

`inc_ricci_v9.cpp` preserves the v8 linear-response and degree-radius tooling, and adds:

- fixed-radius angle-budget diagnostics `B_i`;
- `angle_budget_summary.txt`, `angle_budget_vertex.csv`, and `B_by_degree.csv`;
- fast local edge-flip angle-budget delta evaluation;
- angle-budget preconditioning;
- `angle_mismatch` proposal mode for linear surgery;
- top-K full linear predictor evaluation after cheap local candidate ranking;
- v9 strategy comparison;
- threaded parallel multistart scan.

## Single Angle Precondition

Command:

```bash
./inc_ricci_v9 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --degree_target_mode quantile --angle_precondition 1 --angle_pre_steps 2000 --angle_pre_trials 200 --out angle_pre_v9_8_delta025_ratio5 --progress_every 0
```

| Diagnostic | Initial | Final |
|---|---:|---:|
| `E_B` | `220.455928` | `10.734436` |
| `B_rms` | `2.624738` | `0.579181` |
| `B_max` | `7.137063` | `1.449817` |
| `S_degree` | `75.9699` | `8.2500` |
| degree-radius corr | `-0.030071` | `0.739844` |
| degree match fraction | `0.21875` | `0.46875` |

Linear probe after angle-only preconditioning:

| Diagnostic | Baseline | Angle-only |
|---|---:|---:|
| `S_pred` | `2422.3295` | `331.5171` |
| `rms_K_lin` | `5.300877e-03` | `3.217442e-03` |
| `du_rms` | `0.852836` | `0.280044` |
| `du_max` | `1.993017` | `0.855454` |

Interpretation: angle-budget motifs are a strong signal. They reduce the needed radius response by about `3x`, but not to the strict physical thresholds `du_rms < 0.03`, `du_max < 0.1`.

## Degree + Angle + Linear

Command:

```bash
./inc_ricci_v9 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --pipeline angle_degree_linear --precondition_order degree_then_angle --degree_precondition 1 --degree_pre_steps 1000 --degree_pre_trials 100 --angle_precondition 1 --angle_pre_steps 2000 --angle_pre_trials 200 --linear_surgery 1 --linear_surgery_steps 500 --linear_candidate_pool 200 --linear_full_eval_top 20 --proposal_mode angle_mismatch --out degree_angle_linear_v9_8_delta025_ratio5 --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| accepted linear-surgery flips | 4 |
| `S_total` | `316.1136` |
| `S_pred` | `309.1136` |
| `S_degree` | `7.0000` |
| `E_B` | `8.843935` |
| `B_rms` | `0.525712` |
| `B_max` | `1.682344` |
| `rms_K_lin` | `2.943406e-03` |
| `du_rms` | `0.270346` |
| `du_max` | `0.824261` |
| degree-radius corr | `0.772608` |
| linear predicted good | no |

Interpretation: this exact longer precondition run improved the initial topology strongly, but the final linear surgery got stuck after only 4 accepted moves. It remains outside the small-`du` regime.

## Strategy Compare

Command:

```bash
./inc_ricci_v9 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --v9_strategy_compare 1 --degree_pre_steps 800 --degree_pre_trials 80 --angle_pre_steps 1500 --angle_pre_trials 150 --linear_surgery_steps 300 --linear_candidate_pool 150 --linear_full_eval_top 15 --strategy_out v9_strategy_compare_8.csv --progress_every 0
```

| Strategy | `S_pred` | `du_rms` | `du_max` | `B_rms` | degree corr |
|---|---:|---:|---:|---:|---:|
| none | `2422.3295` | `0.852836` | `1.993017` | `2.624738` | `-0.030071` |
| degree_only | `501.4209` | `0.364391` | `0.987951` | `0.708301` | `0.679700` |
| angle_only | `249.8736` | `0.250544` | `0.705829` | `0.507549` | `0.726847` |
| degree_then_angle | `331.8955` | `0.300437` | `0.775676` | `0.544292` | `0.737946` |
| angle_then_degree | `190.5851` | `0.232473` | `0.546980` | `0.550470` | `0.714759` |
| degree_then_angle_then_linear | `158.0859` | `0.219088` | `0.451381` | `0.513252` | `0.756461` |
| angle_then_degree_then_linear | `159.3856` | `0.220410` | `0.452130` | `0.593112` | `0.702637` |

Best strategy by `S_pred`: `degree_then_angle_then_linear`.

Angle preconditioning beats degree-only preconditioning in this case. Ordering matters mildly: angle-before-degree is better before linear surgery, while both linear-surgery strategies end close.

## Parallel Multistart, Strict Thresholds

Command:

```bash
./inc_ricci_v9 --nx 8 --ny 8 --flips 100 --seed 1000 --radius_seed 101 --same_radii_all_starts 1 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --pipeline angle_degree_linear --precondition_order degree_then_angle --degree_precondition 1 --degree_pre_steps 800 --degree_pre_trials 80 --angle_precondition 1 --angle_pre_steps 1500 --angle_pre_trials 150 --linear_surgery 1 --linear_surgery_steps 300 --linear_candidate_pool 150 --linear_full_eval_top 15 --proposal_mode angle_mismatch --parallel_multistart 1 --threads 16 --multistart_count 100 --multistart_out v9_multistart_8_delta025_ratio5.csv --progress_every 0
```

| Statistic | Value |
|---|---:|
| rows | 100 |
| strict `linear_predicted_good` | 0 / 100 |
| best run | 34 |
| best `S_pred` | `71.5983` |
| best `du_rms` | `0.133565` |
| best `du_max` | `0.329037` |
| min `du_rms` | `0.133565` |
| median `du_rms` | `0.190986` |
| min `du_max` | `0.276890` |
| median `du_max` | `0.445223` |
| corr(`B_rms`, `S_pred`) | `0.5291` |
| corr(degree corr, `S_pred`) | `-0.3205` |

## Parallel Multistart, Relaxed Thresholds

Command used `du_rms_allowed=0.10`, `du_max_allowed=0.25`.

| Statistic | Value |
|---|---:|
| rows | 100 |
| relaxed `linear_predicted_good` | 0 / 100 |
| best run | 27 |
| best `S_pred` | `17.0385` |
| best `du_rms` | `0.137409` |
| best `du_max` | `0.295812` |
| min `du_rms` | `0.137409` |
| median `du_rms` | `0.211416` |
| min `du_max` | `0.295812` |
| median `du_max` | `0.487398` |
| corr(`B_rms`, `S_pred`) | `0.3852` |
| corr(degree corr, `S_pred`) | `-0.4059` |

Even relaxed thresholds were not met. Best `du_max` stayed above `0.25`.

## Larger 12x12 Multistart

Command:

```bash
./inc_ricci_v9 --nx 12 --ny 12 --flips 300 --seed 3000 --radius_seed 202 --same_radii_all_starts 1 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --pipeline angle_degree_linear --precondition_order degree_then_angle --degree_precondition 1 --degree_pre_steps 1500 --degree_pre_trials 120 --angle_precondition 1 --angle_pre_steps 2500 --angle_pre_trials 200 --linear_surgery 1 --linear_surgery_steps 500 --linear_candidate_pool 200 --linear_full_eval_top 20 --proposal_mode angle_mismatch --parallel_multistart 1 --threads 16 --multistart_count 32 --multistart_out v9_multistart_12_delta025_ratio5.csv --progress_every 0
```

| Statistic | Value |
|---|---:|
| rows | 32 |
| `linear_predicted_good` | 0 / 32 |
| best run | 16 |
| best `S_pred` | `99.0925` |
| best `du_rms` | `0.163962` |
| best `du_max` | `0.379362` |
| min `du_rms` | `0.163473` |
| median `du_rms` | `0.202786` |
| min `du_max` | `0.367725` |
| median `du_max` | `0.485621` |
| corr(`B_rms`, `S_pred`) | `0.7208` |
| corr(degree corr, `S_pred`) | `-0.3040` |

The same trend persists at `N=144`: angle/degree topology design improves the predictor but still does not enter the small-perturbation physical regime.

## Conclusion

v9 supports the angle-budget hypothesis partially:

- angle-budget preconditioning is a stronger local topology signal than degree-only matching;
- it reduces `E_B`, `B_rms`, `S_pred`, `du_rms`, and `du_max` by large factors;
- combining degree, angle budget, and linear surgery gives the best strategy comparison result.

But v9 does not yet find a physical INC-like candidate:

- strict `du_rms < 0.03`, `du_max < 0.1` was never reached;
- relaxed `du_rms < 0.10`, `du_max < 0.25` was also not reached in 100 starts;
- 12x12 trends remain negative.

Current interpretation:

```text
angle-budget motifs matter: yes
degree alone is sufficient: no
angle + degree + linear surgery reaches relaxed physical thresholds: no
likely missing ingredient: topology generation from geometric/radical-Delaunay-like local structure, not only abstract edge flips
```

