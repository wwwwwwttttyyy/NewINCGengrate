# INC-Ricci v7 linear-response topology-prediction report

## Build And Tests

Build:

```bash
g++ -O3 -std=c++17 inc_ricci_v7.cpp -o inc_ricci_v7
```

Test command:

```bash
./inc_ricci_v7 --test
```

Result: all v7 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Regular `8x8`, equal radii | PASS | `rmsK0=1.776357e-15`, `du_max=2.79e-22`, `S=8.49e-13` |
| Regular `8x8`, weak lognormal `Delta=0.05` | PASS | `rmsK0=1.699704e-01`, `rmsKlin=3.001208e-04`, `du_max=1.119219e-01` |
| Random `8x8`, 50 flips | PASS | linear summary produced |
| Small linear surgery | PASS | score decreased |
| Linear surgery plus refinement | PASS | refinement ran, no fake success |

## v7 Implementation Summary

`inc_ricci_v7.cpp` keeps the useful v6 infrastructure and adds:

- finite-difference curvature Jacobian `L = dK/du`;
- CG solve for the Tikhonov system `(L^T L + lambda I) du = -L^T K0`;
- linear-response predictor diagnostics;
- predictor score `S_pred`;
- degree-radius correlation and mismatch diagnostics;
- mismatch-guided edge flip proposals;
- `linear_surgery_stats.csv`;
- optional nonlinear refinement after linear surgery;
- topology screening mode;
- linear Pareto scan mode.

## Single Predictor

Command:

```bash
./inc_ricci_v7 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --radius_mode fixed --surgery_metric linear --out linear_predictor_8_f100_delta025_ratio5 --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| `rms_K0` | `2.624738` |
| `max_K0` | `7.137063` |
| `rms_K_lin` | `5.300877e-03` |
| `max_K_lin` | `1.641693e-02` |
| `du_rms` | `0.852836` |
| `du_max` | `1.993017` |
| response ratio | `2.019584e-03` |
| degree-radius corr | `-0.030071` |
| `S_pred` | `2422.041419` |
| singular response | yes |
| linear predicted good | no |

Interpretation: the topology is linearly flattenable only with unphysical radius changes. The linear residual is small, but the required `du` is far beyond the allowed `du_rms=0.03`, `du_max=0.1` thresholds.

## Linear Surgery

Command:

```bash
./inc_ricci_v7 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --radius_mode fixed --surgery_metric linear --proposal_mode mismatch --linear_surgery_steps 500 --linear_surgery_trials_per_step 50 --out linear_surgery_8_f100_delta025_ratio5 --progress_every 0
```

| Diagnostic | Initial | Final |
|---|---:|---:|
| accepted flips | - | 66 |
| `S_pred` | `2422.041419` | `103.915700` |
| `rms_K0` | `2.624738` | `0.711393` |
| `max_K0` | `7.137063` | `2.080223` |
| `rms_K_lin` | `5.300877e-03` | `9.414507e-04` |
| `max_K_lin` | `1.641693e-02` | `1.908924e-03` |
| `du_rms` | `0.852836` | `0.178085` |
| `du_max` | `1.993017` | `0.339674` |
| degree-radius corr | `-0.030071` | `0.718638` |
| singular response | yes | yes |
| linear predicted good | no | no |

Interpretation: mismatch-guided flips are doing the intended topology-selection work. They strongly improve degree-radius correlation and reduce the predicted response cost by about `23x`. But the final topology still requires too much radius motion, so it is not a physical predictor success.

No separate random/highK surgery comparison was run in this required batch.

## Linear Refinement

Command:

```bash
./inc_ricci_v7 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --radius_mode fixed --surgery_metric linear --proposal_mode mismatch --linear_surgery_steps 500 --linear_surgery_trials_per_step 50 --refine_after_linear 1 --refine_du_bound 0.15 --refine_max_iter 50000 --develop 1 --check_crossings 1 --out linear_refine_8_f100_delta025_ratio5 --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| final `max_abs_K` | `0.117290` |
| final `rms_K` | `0.025471` |
| max `|u-u0|` | `0.339674` |
| develop spread | `0.403849` |
| period fit RMS | `0.164848` |
| max edge contact error | `0.165052` |
| nonedge overlap count | 0 |
| nonedge max overlap | 0 |
| edge crossing count | 0 |
| face degenerate count | 0 |
| orientation violations | 0 |
| radius ratio | `3.851547` |
| polydispersity delta | `0.300396` |
| `D_sort` | `5.274282e-03` |
| `D_quantile` | `2.181065e-02` |
| physical candidate flag | no |

Interpretation: nonlinear refinement improves curvature and keeps the packing overlap-free/crossing-free, but it does not reach low curvature or development/contact consistency. It also violates the requested `refine_du_bound=0.15` in the final state, so this is not a physical candidate.

## Topology Screening

Command:

```bash
./inc_ricci_v7 --nx 8 --ny 8 --screen_topologies 1 --screen_count 500 --screen_flips 100 --screen_same_radii 1 --seed 10 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --screen_out topology_screen_8_delta025_ratio5.csv --progress_every 0
```

Result: 500 random topologies screened.

| Quantity | Min | Mean | Max |
|---|---:|---:|---:|
| `S_pred` | `1180.559` | `2274.280` | `4086.573` |
| `rms_K_lin` | `2.859810e-03` | `5.360365e-03` | `1.039924e-02` |
| `du_rms` | `0.569761` | `0.780407` | `0.988547` |
| `du_max` | `1.444683` | `2.099116` | `3.471830` |
| degree-radius corr | `-0.355050` | `0.001669` | `0.401922` |

Good topology fraction:

```text
linear_predicted_good = 0 / 500
```

Best screened topology:

| Diagnostic | Value |
|---|---:|
| sample | 234 |
| seed | 244 |
| `S_pred` | `1180.559371` |
| `rms_K0` | `2.017622` |
| `rms_K_lin` | `2.859810e-03` |
| `du_rms` | `0.569761` |
| `du_max` | `1.497061` |
| degree-radius corr | `-8.375005e-04` |

Correlation between degree-radius correlation and `S_pred`:

```text
corr(degree_radius_corr, S_pred) = -0.368
```

Interpretation: better degree-radius correlation tends to lower predictor score, but random flips almost never produce a topology close to the small-`du` physical regime.

## Linear Pareto

Command:

```bash
./inc_ricci_v7 --nx 8 --ny 8 --flips 100 --seed 20 --radius_seed 1000 --target_dist lognormal --linear_pareto_scan 1 --pareto_delta_list 0.15,0.20,0.25,0.30,0.35,0.40 --pareto_ratio_list 2,3,5,8,12 --pareto_seeds 1,2,3 --linear_surgery_steps 200 --linear_surgery_trials_per_step 30 --pareto_out linear_pareto_v7_8.csv --progress_every 0
```

Result: 90 rows written.

Overall:

| Flag | Count |
|---|---:|
| `linear_predicted_good` | 0 / 90 |
| `refined_physical_candidate` | 0 / 90 |

Best point by `S_pred`:

| Diagnostic | Value |
|---|---:|
| target delta | `0.15` |
| target radius ratio | `3` |
| seed | 2 |
| `S_pred` | `55.403905` |
| `rms_K_lin` | `4.557876e-04` |
| `max_K_lin` | `1.310994e-03` |
| `du_rms` | `0.116095` |
| `du_max` | `0.276010` |
| degree-radius corr | `0.645961` |
| linear predicted good | no |

Best score by target delta:

| delta | best `S_pred` | best `rms_K_lin` | best `du_rms` |
|---:|---:|---:|---:|
| 0.15 | `55.4039` | `4.557876e-04` | `0.116095` |
| 0.20 | `57.9957` | `7.640869e-04` | `0.127047` |
| 0.25 | `61.9167` | `1.066888e-03` | `0.131401` |
| 0.30 | `67.5041` | `7.807282e-04` | `0.135143` |
| 0.35 | `83.2789` | `1.098564e-03` | `0.152806` |
| 0.40 | `70.4396` | `8.221844e-04` | `0.133496` |

Best score by target radius ratio:

| ratio | best `S_pred` | best `rms_K_lin` | best `du_rms` |
|---:|---:|---:|---:|
| 2 | `57.9957` | `7.640869e-04` | `0.127047` |
| 3 | `55.4039` | `4.557876e-04` | `0.116095` |
| 5 | `58.7486` | `5.199163e-04` | `0.124757` |
| 8 | `62.5017` | `7.438817e-04` | `0.130206` |
| 12 | `58.4396` | `6.345008e-04` | `0.125124` |

Interpretation: Pareto improves substantially relative to raw random topology, but all best points still require `du_rms` around `0.11-0.15`, which is several times larger than the default allowed `0.03`. `Delta=0.25`, ratio `5` is not feasible under the current predictor thresholds.

## Conclusion

v7 gives a useful new diagnostic separation:

- The linear system can usually reduce curvature residual strongly.
- The real obstruction is the required radius perturbation, not the residual after solving.
- Mismatch-guided topology surgery works as a topology-selection mechanism: it improved `S_pred`, `du_rms`, `du_max`, and degree-radius correlation dramatically.
- Random topology screening shows good topologies are rare or absent in this flip-only search space for `Delta=0.25`, ratio `5`.
- Nonlinear refinement after linear surgery did not turn the predicted topology into a physical candidate.

Current interpretation:

```text
linear predictor finds small-residual topologies: yes
linear predictor finds small-du physical topologies: no
mismatch-guided topology surgery helps: yes
refinement succeeds after linear search: no
Delta ~0.25, ratio ~5 feasible in current search: no
```

The next step should be more aggressive topology search, not more Ricci-flow tuning: multi-start mismatch surgery, batch flips, annealed acceptance, and explicit degree-radius matching constraints before refinement.
