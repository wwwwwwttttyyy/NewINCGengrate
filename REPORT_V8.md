# INC-Ricci v8 degree-radius topology-design report

## Build And Tests

Build:

```bash
g++ -O3 -std=c++17 inc_ricci_v8.cpp -o inc_ricci_v8
```

Test command:

```bash
./inc_ricci_v8 --test
```

Result: all v8 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Regular `8x8`, equal radii, no degree target | PASS | `rmsK0=1.776357e-15`, `du_max=2.79e-22` |
| Target degree assignment | PASS | `sum_target=384=6N`, `target_corr=0.925326` |
| Small degree precondition | PASS | `S_degree` decreased `75.9699 -> 11.0859` |
| Degree precondition + predictor | PASS | linear summary written |
| Strategy compare | PASS | `strategy_compare_v8_test.csv` produced |

## v8 Implementation Summary

`inc_ricci_v8.cpp` preserves the v7 linear-response core and adds:

- target degree assignment from radius quantiles;
- exact target degree sum constraint `sum d_i*=6N`;
- degree score `S_degree`;
- degree-only topology preconditioning;
- degree-aware linear surgery score `S_total = S_pred + alpha_degree*S_degree`;
- degree-aware mismatch proposal weights;
- degree fraction scan;
- topology strategy comparison;
- target degree output in `target_degrees.csv`.

## Target Degree Assignment

For the main `8x8`, `Delta=0.25`, ratio `5` target:

| Field | Value |
|---|---:|
| target degree histogram | `4:5 5:11 6:32 7:11 8:5` |
| target degree sum | `384 = 6N` |
| target degree-radius corr | `0.925326` |

This confirms that the target assignment is strongly radius-correlated while preserving the torus triangulation degree sum.

## Degree Precondition

Command:

```bash
./inc_ricci_v8 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --degree_target_mode quantile --degree_precondition 1 --degree_pre_steps 1000 --degree_pre_trials 100 --surgery_metric linear --out degree_pre_v8_8_delta025_ratio5 --progress_every 0
```

| Diagnostic | Initial random | After degree precondition |
|---|---:|---:|
| `S_degree` | `75.9699` | `9.5000` |
| `E_deg` | - | `0.7500` |
| max degree error | - | 2 |
| degree match fraction | `0.21875` | `0.43750` |
| degree-radius corr | `-0.03007` | `0.61380` |
| accepted precondition flips | - | 44 |

Linear predictor after degree precondition:

| Diagnostic | Value |
|---|---:|
| `rms_K0` | `0.836101` |
| `max_K0` | `1.995894` |
| `rms_K_lin` | `5.020784e-03` |
| `max_K_lin` | `1.429251e-02` |
| `du_rms` | `0.463191` |
| `du_max` | `1.625949` |
| `S_pred` | `1016.6759` |
| `S_total` | `1026.1759` |
| linear predicted good | no |

Interpretation: degree preconditioning alone dramatically fixes degree-radius mismatch, but it is not enough. The linear response still requires unphysical `du`.

## Degree + Linear Surgery

Command:

```bash
./inc_ricci_v8 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --degree_target_mode quantile --degree_precondition 1 --degree_pre_steps 1000 --degree_pre_trials 100 --surgery_metric linear --proposal_mode mismatch --linear_surgery_use_degree 1 --linear_surgery_steps 500 --linear_surgery_trials_per_step 50 --out degree_linear_v8_8_delta025_ratio5 --progress_every 0
```

| Diagnostic | After degree pre | Final |
|---|---:|---:|
| accepted linear-surgery flips | - | 36 |
| `S_total` | `1026.1759` | `90.9278` |
| `S_pred` | `1016.6759` | `85.2403` |
| `S_degree` | `9.5000` | `5.6875` |
| `rms_K0` | `0.836101` | `0.576490` |
| `rms_K_lin` | `5.020784e-03` | `7.147674e-04` |
| `du_rms` | `0.463191` | `0.154640` |
| `du_max` | `1.625949` | `0.330173` |
| degree match fraction | `0.43750` | `0.53125` |
| degree-radius corr | `0.61380` | `0.74116` |
| linear predicted good | no | no |

Comparison to v7 mismatch-only final:

| Metric | v7 mismatch | v8 degree-aware |
|---|---:|---:|
| `S_pred` | `103.9157` | `85.2403` |
| `du_rms` | `0.178085` | `0.154640` |
| `du_max` | `0.339674` | `0.330173` |
| degree-radius corr | `0.71864` | `0.74116` |

Interpretation: degree-radius design helps. It reduces `S_pred`, `du_rms`, and improves degree-radius correlation relative to v7. But `du_rms` and `du_max` remain above the physical thresholds.

## Refinement

Command with `refine_du_bound=0.15`:

```bash
./inc_ricci_v8 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --degree_target_mode quantile --degree_precondition 1 --degree_pre_steps 1000 --degree_pre_trials 100 --surgery_metric linear --proposal_mode mismatch --linear_surgery_use_degree 1 --linear_surgery_steps 500 --linear_surgery_trials_per_step 50 --refine_after_linear 1 --refine_du_bound 0.15 --develop 1 --check_crossings 1 --out degree_refine_v8_8_delta025_ratio5 --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| final `max_abs_K` | `0.825403` |
| final `rms_K` | `0.315147` |
| max `|u-u0|` | `0.150000` |
| develop spread | `5.617607` |
| period fit RMS | `1.927143` |
| max edge contact error | `2.721848` |
| nonedge overlap count | 3 |
| edge crossing count | 0 |
| radius ratio | `2.848693` |
| polydispersity delta | `0.236232` |
| physical candidate | no |

Command with `refine_du_bound=0.25`:

| Diagnostic | Value |
|---|---:|
| final `max_abs_K` | `0.386927` |
| final `rms_K` | `0.142151` |
| max `|u-u0|` | `0.250000` |
| develop spread | `2.376486` |
| period fit RMS | `0.850498` |
| max edge contact error | `1.208617` |
| nonedge overlap count | 0 |
| edge crossing count | 0 |
| radius ratio | `3.202438` |
| polydispersity delta | `0.244207` |
| physical candidate | no |

Interpretation: relaxing the bound from `0.15` to `0.25` helps curvature and geometry, and removes non-edge overlaps, but still does not reach intrinsic/development/contact tolerances.

## Degree Fraction Scan

Command:

```bash
./inc_ricci_v8 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --degree_fraction_scan 1 --degree_precondition 1 --degree_pre_steps 800 --degree_pre_trials 80 --linear_surgery_use_degree 1 --linear_surgery_steps 300 --linear_surgery_trials_per_step 30 --degree_fraction_out degree_fraction_scan_v8_8.csv --progress_every 0
```

| Profile | `S_total` | `S_pred` | `S_degree` | `du_rms` | `du_max` | degree corr | match frac | predicted? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| weak | `94.8663` | `89.1788` | `5.6875` | `0.158712` | `0.337178` | `0.77690` | `0.53125` | no |
| medium | `69.2390` | `63.8640` | `5.3750` | `0.129360` | `0.285702` | `0.80545` | `0.56250` | no |
| strong | `97.7738` | `90.1488` | `7.6250` | `0.158368` | `0.347208` | `0.77189` | `0.53125` | no |
| very strong | `113.1083` | `105.7958` | `7.3125` | `0.172588` | `0.382170` | `0.82558` | `0.56250` | no |

Best profile: `medium`, which is also the default target degree profile. Stronger degree-radius matching overconstrains the topology and worsens `S_total`.

## Strategy Comparison

Command:

```bash
./inc_ricci_v8 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --compare_topology_strategies 1 --degree_pre_steps 800 --degree_pre_trials 80 --linear_surgery_steps 300 --linear_surgery_trials_per_step 30 --strategy_out strategy_compare_v8_8.csv --progress_every 0
```

| Strategy | `S_total` | `S_pred` | `S_degree` | `du_rms` | `du_max` | degree corr | match frac |
|---|---:|---:|---:|---:|---:|---:|---:|
| random flips only | `2498.2994` | `2422.3295` | `75.9699` | `0.852836` | `1.993017` | `-0.03007` | `0.21875` |
| degree precondition only | `507.4209` | `501.4209` | `6.0000` | `0.364391` | `0.987951` | `0.67970` | `0.50000` |
| linear surgery only | `75.0541` | `68.9916` | `6.0625` | `0.137805` | `0.286250` | `0.81264` | `0.68750` |
| degree precondition + linear surgery | `69.2390` | `63.8640` | `5.3750` | `0.129360` | `0.285702` | `0.80545` | `0.56250` |

Interpretation: degree preconditioning plus degree-aware linear surgery gives the best total score and lowest `du_rms`. Linear surgery alone is close, because mismatch surgery can also discover degree-radius correlation. Degree preconditioning helps, but it is not sufficient by itself.

## Larger 12x12 Run

Command:

```bash
./inc_ricci_v8 --nx 12 --ny 12 --flips 300 --seed 2 --radius_seed 202 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --degree_target_mode quantile --degree_precondition 1 --degree_pre_steps 2000 --degree_pre_trials 150 --surgery_metric linear --proposal_mode mismatch --linear_surgery_use_degree 1 --linear_surgery_steps 800 --linear_surgery_trials_per_step 80 --out degree_linear_v8_12_delta025_ratio5 --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| N | 144 |
| `S_total` | `87.0554` |
| `S_pred` | `81.8887` |
| `S_degree` | `5.1667` |
| `rms_K0` | `0.512959` |
| `rms_K_lin` | `1.492252e-03` |
| `du_rms` | `0.148952` |
| `du_max` | `0.332050` |
| degree match fraction | `0.58333` |
| degree-radius corr | `0.79626` |
| target degree-radius corr | `0.93916` |
| linear predicted good | no |
| runtime | `1190 s` |

Interpretation: the larger system shows the same pattern as `8x8`: degree-radius topology design creates high degree-radius correlation and reduces the linear cost, but the required radius perturbation is still too large.

## Conclusion

v8 supports part of the new hypothesis:

- Degree-radius matching is a real topology-selection signal.
- Degree preconditioning dramatically lowers degree mismatch and turns degree-radius correlation strongly positive.
- Degree-aware linear surgery outperforms random topology and modestly improves over v7 mismatch-only surgery.
- The default medium degree profile is best among the tested profiles.

But v8 does not yet produce a physical INC-like candidate:

- `du_rms` remains around `0.13-0.15`, above the default `0.03` target.
- `du_max` remains around `0.29-0.33`, above the default `0.1` target.
- Nonlinear refinement under bounded perturbations does not reach curvature/development/contact tolerances.

Current interpretation:

```text
degree-radius matching matters: yes
degree preconditioning alone is sufficient: no
degree + linear surgery improves over v7: yes
physical candidate achieved: no
likely missing ingredient: higher-order local motifs or geometric/radical-Delaunay construction beyond degree sequence
```

The next step should bias topology construction toward local geometric motifs, not only degree counts: for example, radical-Delaunay-like initial triangulations from synthetic coordinates, multi-flip moves preserving target degree sequence, and local ring statistics around large/small particles.
