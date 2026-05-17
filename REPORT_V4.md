# INC-Ricci v4 packing-diagnostic report

## Build

```bash
g++ -O3 -std=c++17 inc_ricci_v4.cpp -o inc_ricci_v4
```

v4 preserves the v3 fast Ricci/develop path and adds disk-packing diagnostics:

- non-edge overlap checking with periodic minimum image;
- edge contact error checking;
- face signed-area diagnostics;
- optional periodic 3x3 edge crossing checks;
- radius distribution statistics;
- bounds scan over `u in [-B,B]`;
- small-system greedy topology surgery.

## Tests

Command:

```bash
./inc_ricci_v4 --test
```

Result: all tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Regular `8x8`, `u=0`, overlap and crossings | PASS | `maxK=1.776357e-15`, `nonedge_overlap=0`, `crossings=0`, `physical=1` |
| Regular `8x8`, random `u`, unbounded Ricci | PASS | `finalK=8.527936e-10`, `nonedge_overlap=0`, `ratio=1` |
| Random `8x8`, 50 flips, unbounded | PASS | `finalK=7.851120e-10`, `nonedge_overlap=0`, `valid_disk=0` |
| Random `8x8`, 50 flips, bounded `B=2` | PASS | `finalK=7.851120e-10`, no fake success required |
| Small bounds scan | PASS | `bounds_scan_test.csv` produced |

## Bounds Scan

Command:

```bash
./inc_ricci_v4 --nx 8 --ny 8 --flips 100 --seed 1 --method ricci --max_iter 100000 --tol 1e-9 --develop 1 --scan_bounds 1 --bounds_list 0.5,1.0,1.5,2.0,2.5,3.0,4.0 --include_unbounded 1 --bounds_out bounds_scan_8_f100.csv --progress_every 0 --scan_write_details 0
```

Key findings:

- Smallest bounded `B` reaching intrinsic convergence: `B=2.0`.
- Smallest bounded `B` reaching developed consistency: `B=2.0`.
- Smallest bounded `B` with no non-edge overlap: `B=2.0`.
- At `B=2.0`, the result is not physically reasonable by default radius thresholds.

| B | final `max_abs_K` | radius ratio | delta | nonedge overlaps | max overlap | physical radius? |
|---:|---:|---:|---:|---:|---:|---|
| unbounded | `9.238463e-10` | `22.2627` | `0.6540` | 0 | 0 | no |
| 0.5 | `3.329798e+00` | `2.6040` | `0.2139` | 218 | `1.8950` | yes |
| 1.0 | `1.164787e+00` | `7.3815` | `0.5248` | 25 | `1.5917` | no |
| 1.5 | `3.929845e-01` | `14.2109` | `0.6313` | 2 | `0.3971` | no |
| 2.0 | `9.238463e-10` | `22.2627` | `0.6540` | 0 | 0 | no |

Interpretation: tight bounds keep radii closer to physical scale but block curvature closure and create overlaps. Bounds large enough to close curvature allow broad radii.

## Flip Scan

Command:

```bash
./inc_ricci_v4 --nx 8 --ny 8 --method ricci --max_iter 100000 --tol 1e-9 --develop 1 --scan_flips 1 --scan_list 0,1,5,10,50,100,500 --scan_seeds 1,2,3 --bounded 0 --check_crossings 0 --scan_out scan_v4_unbounded_8.csv --progress_every 0 --scan_write_details 0
```

Aggregated results:

| flips | max nonedge overlap | avg radius ratio | max radius ratio | valid disk count | physical radius count |
|---:|---:|---:|---:|---:|---:|
| 0 | 0 | `1.00` | `1.00` | 3 | 3 |
| 1 | 0 | `1.66` | `1.66` | 2 | 3 |
| 5 | 0 | `2.79` | `3.30` | 0 | 3 |
| 10 | 0 | `3.60` | `4.14` | 0 | 3 |
| 50 | 0 | `19.62` | `23.13` | 0 | 0 |
| 100 | 0 | `42.54` | `64.04` | 0 | 0 |
| 500 | 0 | `771.08` | `1863.07` | 0 | 0 |

Interpretation: in these `8x8` unbounded runs, non-edge overlaps were absent, but radius ratio grows rapidly with flip count. The algorithm increasingly produces mathematical circle-packing metrics with physically unreasonable polydispersity.

## Medium Detailed Run

Command:

```bash
./inc_ricci_v4 --nx 16 --ny 16 --flips 2000 --seed 1 --max_iter 200000 --tol 1e-10 --method ricci --bounded 0 --develop 1 --check_crossings 1 --out run_16_v4_unbounded_detailed --progress_every 1000
```

| Diagnostic | Value |
|---|---:|
| final `max_abs_K` | `9.90412196699708147e-11` |
| radius ratio | `1.91425226572332349e+03` |
| polydispersity delta | `1.18505605526331714e+00` |
| develop spread | `1.60766112012541555e-08` |
| period fit RMS | `2.12822663316447427e-08` |
| edge contact max error | `2.25049756608086682e-08` |
| nonedge overlap count | 0 |
| nonedge max overlap | 0 |
| edge crossing count | 0 |
| face bad area count | 242 |
| valid disk packing candidate | no |
| strong INC candidate | no |
| physical INC-like candidate | no |

Interpretation: curvature closes, the developed metric is consistent, edge contacts are accurate, and no non-edge overlaps were detected. But the signed-area diagnostic finds many flipped/folded faces, and the radius distribution is far too broad. This is not a physical INC-like packing.

## Surgery Test

Command:

```bash
./inc_ricci_v4 --nx 8 --ny 8 --flips 100 --seed 1 --max_iter 50000 --tol 1e-8 --method ricci --bounded 1 --u_min -2 --u_max 2 --develop 1 --check_crossings 0 --surgery 1 --surgery_steps 100 --surgery_trials_per_step 10 --surgery_metric mixed --surgery_ricci_iter 10000 --out surgery_8_test --progress_every 0
```

Result:

- Accepted surgery moves: 63.
- Score decreased from `4.12053855736337393e+00` to `7.64597781991114005e-05`.
- Radius ratio improved from about `22.79` at step 1 to `3.10`.
- Polydispersity delta improved to `0.2343`.
- Non-edge overlaps stayed at 0.
- Curvature remained low: final `max_abs_K = 3.93988419666868595e-09`.
- The run still failed `valid_disk_packing_candidate` because `face_bad_area_count = 83`.

Interpretation: surgery is worth developing further. It significantly improves radius distribution under bounded radii while preserving low curvature and overlap-free status, but the current score does not penalize signed-area/folding strongly enough.

## Overall Conclusion

v4 changes the interpretation:

- Intrinsic Ricci closure alone is not enough.
- Many unbounded random-topology endpoints are valid as mathematical intrinsic circle-packing metrics, but not physical INC-like packings because radius ratios become huge.
- Bounded radii expose a tradeoff: tight physical bounds block convergence or introduce overlaps, while looser bounds recover convergence at the cost of broad radii.
- Topology surgery appears to be the most promising next step because it can reduce radius ratio and polydispersity substantially.

Current status:

```text
intrinsic endpoint: achieved in many runs
overlap-free disk metric: sometimes achieved
strong INC candidate: not achieved in detailed randomized medium run
physical INC-like candidate: not achieved in randomized production run
```

The next useful v5 step is to include face orientation/folding and radius distribution directly in the surgery score, then rerun bounded small and medium systems.
