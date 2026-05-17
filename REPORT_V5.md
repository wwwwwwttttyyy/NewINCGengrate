# INC-Ricci v5 bounded-surgery report

## Build And Tests

Build:

```bash
g++ -O3 -std=c++17 inc_ricci_v5.cpp -o inc_ricci_v5
```

Test command:

```bash
./inc_ricci_v5 --test
```

Result: all v5 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Regular `8x8`, `u=0` | PASS | `maxK=1.776357e-15`, overlap-free, no crossings, physical flag true |
| Regular `8x8`, random `u`, unbounded | PASS | `finalK=8.527936e-10`, radius ratio `1` |
| Random `8x8`, 50 flips, unbounded | PASS | `finalK=7.851120e-10`, program reports candidate flags |
| Random `8x8`, 50 flips, bounded `[-1,1]` | PASS | bounded flow stalls with `finalK=1.179827e+00`, as expected |
| Small surgery | PASS | score decreased from `7.080680e+01` to `5.360118e-01` |
| Small multistart | PASS | multistart CSV produced |

## What Changed In v5

`inc_ricci_v5.cpp` keeps the v4 fast Ricci/develop/diagnostic path and adds:

- separated face diagnostics: `face_degenerate_count`, `face_negative_signed_count`, and `orientation_neighbor_violation_count`;
- `--require_consistent_orientation` and `--check_orientation_neighbors`;
- v5 surgery score with explicit curvature, development, overlap, crossing, degeneracy, orientation, radius-ratio, and polydispersity terms;
- bounded surgery that rejects Ricci steps outside `[u_min,u_max]` without clamping;
- multistart surgery output;
- disorder diagnostics: degree histogram, degree-6 fraction, edge-length CV, radius CV, and a lightweight noncrystalline topology heuristic.

Candidate flags remain separate:

```text
valid_disk_packing_candidate
strong_INC_candidate
physical_INC_like_candidate
```

## Bounded Surgery, B=1

Command:

```bash
./inc_ricci_v5 --nx 8 --ny 8 --flips 100 --seed 1 --method ricci --max_iter 50000 --tol 1e-8 --bounded 1 --u_min -1.0 --u_max 1.0 --develop 1 --check_crossings 1 --target_max_ratio 5 --target_delta 0.25 --surgery 1 --surgery_steps 300 --surgery_trials_per_step 30 --surgery_metric v5 --surgery_ricci_iter 10000 --out surgery_v5_8_f100_B1 --progress_every 0
```

| Diagnostic | Initial | Final |
|---|---:|---:|
| score | `2.279364e+02` | `2.173913e+01` |
| accepted surgery moves | - | 44 |
| `max_abs_K` | `1.541404e+00` | `1.543593e+00` |
| develop spread | `1.333028e+01` | `3.234938e+00` |
| period fit RMS | `3.432387e+00` | `1.177462e+00` |
| nonedge max overlap | `1.645748e+00` | `0` |
| nonedge overlap count | 24 | 0 |
| edge crossing count | 6 | 0 |
| face degenerate count | 0 | 0 |
| orientation neighbor violations | 7 | 0 |
| radius ratio | `7.381509` | `7.141672` |
| polydispersity delta | `0.524802` | `0.508569` |
| physical INC-like candidate | no | no |

Interpretation: surgery removed overlaps, crossings, and neighbor-orientation violations, but bounded `B=1` did not reach low curvature, development consistency, contact validity, or physical radius statistics.

## Bounded Surgery, B=1.5

Command:

```bash
./inc_ricci_v5 --nx 8 --ny 8 --flips 100 --seed 1 --method ricci --max_iter 50000 --tol 1e-8 --bounded 1 --u_min -1.5 --u_max 1.5 --develop 1 --check_crossings 1 --target_max_ratio 8 --target_delta 0.35 --surgery 1 --surgery_steps 300 --surgery_trials_per_step 30 --surgery_metric v5 --surgery_ricci_iter 10000 --out surgery_v5_8_f100_B15 --progress_every 0
```

| Diagnostic | Initial | Final |
|---|---:|---:|
| score | `6.345399e+01` | `2.628034e+01` |
| accepted surgery moves | - | 6 |
| `max_abs_K` | `1.781045e+00` | `1.678105e+00` |
| develop spread | `4.359440e+00` | `4.713166e+00` |
| period fit RMS | `1.480421e+00` | `2.167432e+00` |
| nonedge max overlap | `2.733264e-01` | `0` |
| nonedge overlap count | 2 | 0 |
| edge crossing count | 0 | 0 |
| face degenerate count | 0 | 0 |
| orientation neighbor violations | 1 | 1 |
| radius ratio | `14.210923` | `14.210923` |
| polydispersity delta | `0.631288` | `0.631288` |
| physical INC-like candidate | no | no |

Interpretation: `B=1.5` also failed the physical target. It removed non-edge overlap but kept large curvature/development residuals and poor radius statistics.

## Multistart Surgery, B=1

Command:

```bash
./inc_ricci_v5 --nx 8 --ny 8 --flips 100 --seed 10 --method ricci --max_iter 50000 --tol 1e-8 --bounded 1 --u_min -1.0 --u_max 1.0 --develop 1 --check_crossings 0 --target_max_ratio 5 --target_delta 0.25 --multistart_surgery 1 --multistart_count 10 --surgery_steps 200 --surgery_trials_per_step 20 --surgery_metric v5 --surgery_ricci_iter 10000 --multistart_out multistart_v5_B1.csv --progress_every 0
```

Success rates:

| Flag | Count |
|---|---:|
| `valid_disk_packing_candidate` | 5 / 10 |
| `strong_INC_candidate` | 5 / 10 |
| `physical_INC_like_candidate` | 0 / 10 |

Distribution:

| Quantity | Min | Mean | Max |
|---|---:|---:|---:|
| radius ratio | `4.849198` | `6.009120` | `7.305337` |
| `max_abs_K` | `6.252822e-09` | `7.176454e-01` | `2.364790` |

Best rows:

| run | seed | score | `max_abs_K` | radius ratio | delta | physical? |
|---:|---:|---:|---:|---:|---:|---|
| 0 | 10 | `4.121494e-01` | `6.252822e-09` | `5.280208` | `0.403227` | no |
| 6 | 16 | `4.595062e-01` | `7.471278e-09` | `5.129545` | `0.412206` | no |
| 9 | 19 | `4.608483e-01` | `6.354710e-09` | `4.849198` | `0.426642` | no |

Interpretation: multistart bounded surgery sometimes finds coordinate-level strong candidates under the geometric checks used here, but none satisfy the physical radius distribution target. The limiting quantity is mostly polydispersity, and sometimes radius ratio.

## Unbounded Baseline

Command:

```bash
./inc_ricci_v5 --nx 8 --ny 8 --flips 100 --seed 1 --method ricci --max_iter 50000 --tol 1e-9 --bounded 0 --develop 1 --check_crossings 1 --target_max_ratio 5 --target_delta 0.25 --out baseline_unbounded_8_f100 --progress_every 0
```

| Diagnostic | Value |
|---|---:|
| final `max_abs_K` | `9.238463e-10` |
| develop spread | `1.538960e-08` |
| period fit RMS | `1.834050e-08` |
| theta mean abs | `3.222972e-16` |
| max edge contact error | `1.811472e-08` |
| nonedge overlap count | 0 |
| nonedge max overlap | 0 |
| edge crossing count | 0 |
| face degenerate count | 0 |
| orientation neighbor violations | 0 |
| radius ratio | `22.262693` |
| polydispersity delta | `0.654016` |
| valid disk packing candidate | yes |
| strong INC candidate | yes |
| physical INC-like candidate | no |

Interpretation: the unbounded baseline is a valid coordinate-level mathematical disk-packing candidate under the implemented diagnostics, but not physical INC-like because the radius distribution is too broad.

## Conclusion

v5 confirms the core distinction:

- intrinsic Ricci closure is not the bottleneck;
- coordinate-level strong candidates can appear in unbounded or some bounded multistart cases;
- physical INC-like candidates were not found in these required runs because bounded radii either block curvature/development closure or leave polydispersity too high;
- topology surgery is useful because it removes overlaps/crossings and improves scores, but the current v5 score still needs stronger pressure toward low polydispersity without sacrificing curvature and development consistency.

Current status:

```text
intrinsic endpoint: achieved in unbounded baseline and some multistart bounded runs
coordinate-level strong candidate: achieved in unbounded baseline and 5/10 multistart bounded runs
physical INC-like candidate: not achieved in the required v5 experiments
```

Recommended next step: tune surgery around radius statistics directly, especially `wRad`, `wDelta`, and topology proposal size, and run larger multistart searches at bounds near `B=1.0` to `B=1.5`.
