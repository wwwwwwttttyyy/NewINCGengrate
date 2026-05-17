# INC-Ricci v3 performance report

## Build

```bash
g++ -O3 -std=c++17 inc_ricci_v3.cpp -o inc_ricci_v3
```

v3 focuses on the fast production path:

- Ricci curvature computation is O(F) per iteration.
- Topology, edge adjacency, face adjacency, incident faces, and canonical periodic shifts are built after topology generation, not inside Ricci iterations.
- Temporary edge maps use packed `uint64_t` keys.
- Scan mode writes one CSV row per run and uses up to 4 worker threads by default.
- Development diagnostics are O(F): BFS over faces once, spread grouping by lifted vertex key, and period vector fitting with one reference occurrence per base vertex.

## Tests

Command:

```bash
./inc_ricci_v3 --test
```

Result: PASS.

| Test | Result | Key diagnostics |
|---|---|---|
| Regular `8x8` flat develop | PASS | `maxK=1.776357e-15`, `theta_mean=1.572816e-16`, `period_fit=2.348255e-15` |
| Random `8x8` Ricci/develop | PASS | `flips=172`, `finalK=9.241204e-09`, `spread=1.290595e-07` |
| Bounded random `8x8` Ricci | PASS | `finalK=2.141768e-01`, `boundary_rejections=46` |

## Fast Scan

Command:

```bash
./inc_ricci_v3 --nx 8 --ny 8 --method ricci --max_iter 100000 --tol 1e-9 --develop 1 --scan_flips 1 --scan_list 0,1,5,10,50,100,500 --scan_seeds 1,2,3 --bounded 0 --scan_out scan_unbounded_8.csv --progress_every 0 --scan_write_details 0
```

Result:

- rows: 21
- threads: 4
- total wall time: `2.789150e-02 s`
- output: `scan_unbounded_8.csv`

All scan rows completed with `status=ok`. The slowest rows were the `flips=500` cases at about `1.4e-2` to `1.6e-2 s` each.

## Medium Benchmark

Command:

```bash
./inc_ricci_v3 --nx 16 --ny 16 --flips 2000 --seed 1 --max_iter 200000 --tol 1e-10 --method ricci --bounded 0 --develop 1 --embed 0 --out run_16_unbounded_develop --progress_every 1000
```

| Diagnostic | Value |
|---|---:|
| N, E, F | `256, 768, 512` |
| flips requested / accepted | `2000 / 1600` |
| Ricci iterations | `564` |
| final `E_K` | `2.10831705536706334e-19` |
| final `max_abs_K` | `9.90412196699708147e-11` |
| radius ratio | `1.91425226572332349e+03` |
| develop global vertex max spread | `1.60766112012541555e-08` |
| period fit RMS | `2.12822663316447427e-08` |
| develop theta mean abs | `1.92882276990602453e-15` |
| topology time | `1.320388e-01 s` |
| Ricci time | `1.754000e-02 s` |
| develop time | `3.161000e-04 s` |
| total time | `1.499112e-01 s` |

## Large Benchmark

Command:

```bash
./inc_ricci_v3 --nx 32 --ny 32 --flips 5000 --seed 1 --max_iter 300000 --tol 1e-9 --method ricci --bounded 0 --develop 1 --embed 0 --out run_32_unbounded_develop --progress_every 2000
```

| Diagnostic | Value |
|---|---:|
| N, E, F | `1024, 3072, 2048` |
| flips requested / accepted | `5000 / 4247` |
| Ricci iterations | `1156` |
| final `E_K` | `9.86208631776850548e-17` |
| final `max_abs_K` | `9.93273907567981951e-10` |
| radius ratio | `5.01232616450209150e+02` |
| develop global vertex max spread | `3.17275972699945160e-06` |
| period fit RMS | `1.95031158299349645e-06` |
| develop theta mean abs | `1.61295966636584397e-15` |
| topology time | `1.475869e+00 s` |
| Ricci time | `1.449860e-01 s` |
| develop time | `7.651000e-04 s` |
| total time | `1.621671e+00 s` |

## Bottleneck

The dominant cost is topology generation by random edge flips.

For the `32x32` large run:

- topology/flips: about `1.476 s`
- Ricci flow: about `0.145 s`
- linear development diagnostics: about `0.0008 s`

The Ricci core is no longer the bottleneck. The remaining bottleneck is that v3 rebuilds the temporary packed edge map during flip attempts. This is acceptable for the requested `nx<=32` benchmarks, but larger scans should move to an incremental edge map or batched flip strategy.

## Interpretation

The v3 benchmarks show fast intrinsic convergence and fast linear development diagnostics. The reported `develop_theta_mean_abs` is near machine precision because each developed face is constructed from the intrinsic tangent triangle lengths.

This is not the same as a full coordinate-level packing solver. The development spread and period-fit diagnostics are lightweight tests of global consistency; they do not replace the stronger stress embedding / coordinate packing test from v2.
