# INC-Ricci v6 distribution-preserving topology-selection report

## Build And Tests

Build:

```bash
g++ -O3 -std=c++17 inc_ricci_v6.cpp -o inc_ricci_v6
```

Test command:

```bash
./inc_ricci_v6 --test
```

Result: all v6 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Regular `8x8`, fixed equal radii | PASS | `maxK=1.776357e-15`, `physical_v6=1`, crystalline heuristic false |
| Random `8x8`, 50 flips, fixed lognormal target | PASS | program runs, `rmsK=1.998504` |
| Random `8x8`, weak mode | PASS | distribution penalties computed, `D_sort=1.125225e-02` |
| Small fixed surgery | PASS | score decreased `1.595626e+02 -> 9.806527e+01` |
| Small Pareto | PASS | `pareto_v6_test.csv` produced with 8 full-factorial rows |

## v6 Implementation Summary

`inc_ricci_v6.cpp` preserves the v5 intrinsic Ricci/develop/packing diagnostics and adds:

- `--radius_mode free|fixed|weak`;
- target radius distributions: `lognormal`, `uniform_u`, `powerlaw_like`;
- `target_radii.csv` and `distribution_stats.txt`;
- fixed-radius curvature evaluation with `fixed_intrinsic_good`;
- weak distribution-preserving Ricci-like relaxation using `E_reg`;
- v6 surgery score with curvature, RMS curvature, development, overlap, crossing, orientation, `D_sort`, delta, and radius-ratio terms;
- proposal modes `random`, `highK`, and `overlap`;
- Pareto scan mode over target delta, radius ratio, radius mode, and seeds;
- stricter `physical_INC_like_candidate_v6`.

## Fixed-Radius Surgery

Command:

```bash
./inc_ricci_v6 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --radius_mode fixed --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --method ricci --max_iter 1 --develop 1 --check_crossings 1 --target_max_ratio 5 --surgery 1 --surgery_steps 500 --surgery_trials_per_step 40 --surgery_metric v6 --proposal_mode highK --out fixed_surgery_v6_8_f100_delta025_ratio5 --progress_every 0
```

| Diagnostic | Initial | Final |
|---|---:|---:|
| score | `4.515690e+04` | `4.176400e+03` |
| accepted surgery moves | - | 47 |
| `max_abs_K` | `7.137063` | `5.135574` |
| `rms_K` | - | `2.276056` |
| develop spread | `12.892260` | `16.048328` |
| period fit RMS | `4.664486` | `5.315164` |
| nonedge overlap count | 791 | 159 |
| nonedge max overlap | `2.467865` | `2.172439` |
| edge crossing count | 3125 | 400 |
| face degenerate count | 0 | 0 |
| orientation neighbor violations | 31 | 34 |
| radius ratio | `2.814205` | `2.814205` |
| polydispersity delta | `0.244376` | `0.244376` |
| `D_sort` | near 0 | near 0 |
| `fixed_intrinsic_good` | no | no |
| `strong_INC_candidate` | no | no |
| `physical_INC_like_candidate_v6` | no | no |

Interpretation: fixed-radius surgery preserves a physical radius distribution, but cannot make this topology legal or low-curvature. It reduces crossings and overlaps substantially, but the final state still has large curvature, inconsistent development, non-edge overlaps, and crossings.

## Weak Distribution-Preserving Surgery

Command:

```bash
./inc_ricci_v6 --nx 8 --ny 8 --flips 100 --seed 1 --radius_seed 101 --radius_mode weak --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --method ricci --max_iter 50000 --tol 1e-8 --develop 1 --check_crossings 1 --target_max_ratio 5 --bounded 1 --u_min -1.0 --u_max 1.0 --lambda_sort 200 --lambda_delta 100 --lambda_ratio 100 --surgery 1 --surgery_steps 500 --surgery_trials_per_step 40 --surgery_metric v6 --proposal_mode highK --surgery_ricci_iter 10000 --out weak_surgery_v6_8_f100_delta025_ratio5 --progress_every 0
```

| Diagnostic | Initial | Final |
|---|---:|---:|
| score | `4.149395e+03` | `1.207752e+02` |
| accepted surgery moves | - | 37 |
| `max_abs_K` | `1.842153` | `3.059368` |
| `rms_K` | - | `0.928234` |
| develop spread | `16.814025` | `5.214902` |
| period fit RMS | `5.482233` | `2.400240` |
| nonedge overlap count | 100 | 12 |
| nonedge max overlap | `3.059504` | `0.680193` |
| edge crossing count | 189 | 0 |
| face degenerate count | 0 | 0 |
| orientation neighbor violations | 23 | 2 |
| radius ratio | `5.719148` | `7.166811` |
| polydispersity delta | `0.423537` | `0.496480` |
| `D_sort` | `0.040200` | `0.073516` |
| `D_quantile` | `0.157257` | `0.307260` |
| `D_delta` | `0.030115` | `0.060752` |
| `D_ratio` | `0.018058` | `0.129616` |
| `weak_distribution_preserved` | no | no |
| `strong_INC_candidate` | no | no |
| `physical_INC_like_candidate_v6` | no | no |

Interpretation: weak mode improves geometry more than fixed mode and removes crossings, but it does so by drifting away from the target distribution. It still has large curvature, development residuals, overlaps, and orientation violations.

## Pareto Scan

Command:

```bash
./inc_ricci_v6 --nx 8 --ny 8 --flips 100 --seed 10 --radius_seed 1000 --target_dist lognormal --develop 1 --check_crossings 0 --surgery 1 --surgery_steps 200 --surgery_trials_per_step 20 --surgery_metric v6 --proposal_mode highK --pareto_scan 1 --pareto_delta_list 0.20,0.25,0.30,0.35,0.40 --pareto_ratio_list 3,5,8,12 --pareto_modes fixed,weak --pareto_seeds 1,2,3 --pareto_out pareto_v6_8_f100.csv --progress_every 0
```

Result: 120 rows written. No successful v6 physical candidate was found.

Overall success counts:

| Flag | Count |
|---|---:|
| `fixed_intrinsic_good` | 0 / 120 |
| `valid_disk_packing_candidate` | 0 / 120 |
| `strong_INC_candidate` | 0 / 120 |
| `physical_INC_like_candidate_v6` | 0 / 120 |

Success by target delta:

| target delta | runs | fixed good | valid disk | strong | physical v6 | best score |
|---:|---:|---:|---:|---:|---:|---:|
| 0.20 | 24 | 0 | 0 | 0 | 0 | `63.3054` |
| 0.25 | 24 | 0 | 0 | 0 | 0 | `70.0287` |
| 0.30 | 24 | 0 | 0 | 0 | 0 | `19.2979` |
| 0.35 | 24 | 0 | 0 | 0 | 0 | `21.8968` |
| 0.40 | 24 | 0 | 0 | 0 | 0 | `19.2099` |

Success by target radius ratio:

| target ratio | runs | fixed good | valid disk | strong | physical v6 | best score |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 30 | 0 | 0 | 0 | 0 | `60.4751` |
| 5 | 30 | 0 | 0 | 0 | 0 | `60.2171` |
| 8 | 30 | 0 | 0 | 0 | 0 | `19.2099` |
| 12 | 30 | 0 | 0 | 0 | 0 | `19.2979` |

Best tradeoff row by v6 score:

| Field | Value |
|---|---:|
| mode | `weak` |
| target delta | `0.40` |
| target radius ratio | `8` |
| seed | `2` |
| score | `19.209851` |
| `max_abs_K` | `1.201274e-07` |
| `rms_K` | `2.974397e-08` |
| develop spread | `6.425536e-07` |
| period fit RMS | `3.411734e-07` |
| nonedge overlap count | 0 |
| nonedge max overlap | 0 |
| edge crossing count | 0 |
| orientation violations | 0 |
| radius ratio | `8.163241` |
| polydispersity delta | `0.533754` |
| `D_sort` | `0.042835` |
| `D_quantile` | `0.187942` |
| `physical_INC_like_candidate_v6` | no |

Interpretation: the best Pareto point is close geometrically, but it still fails the strict v6 classification because curvature is above `tol=1e-8` and the radius distribution is not preserved. The frontier improves when the target distribution is broader, especially around `target_delta >= 0.30` and `target_radius_ratio >= 8`, but not enough to count as physical.

## Free Baseline

Command:

```bash
./inc_ricci_v6 --nx 8 --ny 8 --flips 100 --seed 1 --radius_mode free --method ricci --max_iter 50000 --tol 1e-9 --develop 1 --check_crossings 1 --target_max_ratio 5 --target_delta 0.25 --out free_baseline_v6_8_f100 --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| `max_abs_K` | `9.238463e-10` |
| `rms_K` | `4.142068e-10` |
| develop spread | `1.538960e-08` |
| period fit RMS | `1.834050e-08` |
| theta mean abs | `3.222972e-16` |
| nonedge overlap count | 0 |
| edge crossing count | 0 |
| radius ratio | `22.262693` |
| polydispersity delta | `0.654016` |
| valid disk packing candidate | yes |
| strong INC candidate | yes |
| physical INC-like candidate | no |
| physical INC-like candidate v6 | no |

Interpretation: free mode still produces a coordinate-level strong mathematical circle-packing candidate, but only by broadening the radius distribution far beyond the target.

## Conclusion

v6 shifts the test from free radius closure to distribution-preserving topology selection. The result is negative for the required experiments:

- fixed radii preserve the target distribution but do not achieve low curvature or legal geometry;
- weak radii improve geometry but drift away from the target distribution;
- Pareto scan found no fixed, strong, or physical v6 candidates;
- free mode remains successful only as a mathematical circle-packing construction with broad radii.

Current interpretation:

```text
fixed mode succeeds: no
weak mode succeeds: no
free mode succeeds: yes, but not physically
topology-selection hypothesis: not supported by this v6 search at delta=0.25, ratio=5
best observed tradeoff: weak mode, broader target distribution near delta=0.40 and ratio=8
```

The next useful step is not to add dynamics, but to improve topology search: larger multi-start Pareto scans, stronger local proposal mechanisms around curvature/overlap defects, and possibly accepting multi-flip topology moves that can cross score barriers.
