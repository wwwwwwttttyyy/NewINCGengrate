# INC-Ricci v2 coordinate-level diagnostic report

## Build

```bash
g++ -O3 -std=c++17 inc_ricci_v2.cpp -o inc_ricci_v2
```

`inc_ricci_v2.cpp` preserves the v1 intrinsic Ricci-flow core and adds:

- bounded log-radius Ricci steps with rejection instead of clamping;
- lifted periodic face corners and canonical edge shifts;
- edge-shift validation after topology generation and flips;
- coordinate embedding stress optimization;
- measured coordinate-level theta residuals;
- strong-INC candidate classification requiring intrinsic closure, edge-length embedding, and small coordinate theta residual.

## Tests

Command:

```bash
./inc_ricci_v2 --test
```

Result: all tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Regular `8x8`, `u_i=0`, embedding | PASS | `max_abs_K=1.776357e-15`, `max_edge_res=2.220446e-16`, `theta_mean_abs=1.572816e-16` |
| Random `8x8`, unbounded Ricci | PASS | `initial_max_abs_K=8.852255e+00`, `final_max_abs_K=8.686049e-09` |
| Random `8x8`, bounded Ricci `[-2,2]` | PASS | `final_max_abs_K=4.578783e-01`, `boundary_rejections=45` |
| Embedding diagnostic of bounded test result | PASS | `E_embed=1.856780e+02`, `max_edge_res=5.526244e+00`, `theta_mean_abs=5.018625e-01` |

## Production Runs

### `run_16_unbounded_embed`

Command:

```bash
./inc_ricci_v2 --nx 16 --ny 16 --flips 2000 --seed 1 --max_iter 200000 --tol 1e-10 --method ricci --bounded 0 --embed 1 --embed_iter 200000 --embed_tol 1e-8 --out run_16_unbounded_embed
```

| Diagnostic | Value |
|---|---:|
| final `E_K` | `1.55527529116044540e-19` |
| final `max_abs_K` | `9.80160397290319452e-11` |
| radius ratio | `3.20105755987180601e+02` |
| boundary rejections | `0` |
| final `max_edge_abs_residual` | `2.35745304171754491e+01` |
| final `theta_mean_abs` | `5.08499999429680272e-01` |
| final `theta_max_abs` | `2.84050372565739728e+00` |
| strong-INC candidate | no |

Interpretation: intrinsic Ricci closure succeeds, but the coordinate embedding fails badly. This is only an intrinsic endpoint.

### `run_16_bounded_embed`

Command:

```bash
./inc_ricci_v2 --nx 16 --ny 16 --flips 2000 --seed 1 --max_iter 200000 --tol 1e-8 --method ricci --bounded 1 --u_min -2 --u_max 2 --embed 1 --embed_iter 200000 --embed_tol 1e-8 --out run_16_bounded_embed
```

| Diagnostic | Value |
|---|---:|
| final `E_K` | `3.75554036877481607e+01` |
| final `max_abs_K` | `1.33961338489385651e+00` |
| radius ratio | `2.74352000817989143e+01` |
| boundary rejections | `44` |
| final `max_edge_abs_residual` | `1.15009112914223337e+01` |
| final `theta_mean_abs` | `4.91245918475712584e-01` |
| final `theta_max_abs` | `2.64276041808047468e+00` |
| strong-INC candidate | no |

Interpretation: bounded radii prevent extreme radius growth, but the `[-2,2]` box blocks Ricci convergence for this topology and seed. Coordinate embedding also fails.

### `run_12_bounded_embed`

Command:

```bash
./inc_ricci_v2 --nx 12 --ny 12 --flips 1000 --seed 2 --max_iter 200000 --tol 1e-8 --method ricci --bounded 1 --u_min -2.5 --u_max 2.5 --embed 1 --embed_iter 200000 --embed_tol 1e-8 --out run_12_bounded_embed
```

| Diagnostic | Value |
|---|---:|
| final `E_K` | `6.78501497919811736e-01` |
| final `max_abs_K` | `2.79462111059942675e-01` |
| radius ratio | `1.16768889363680401e+02` |
| boundary rejections | `46` |
| final `max_edge_abs_residual` | `1.11915127162903154e+01` |
| final `theta_mean_abs` | `4.99287417465770877e-01` |
| final `theta_max_abs` | `2.42072428142736440e+00` |
| strong-INC candidate | no |

Interpretation: the wider `[-2.5,2.5]` box gets closer intrinsically than the `16x16` bounded case, but still stalls before the requested curvature tolerance. Coordinate embedding does not realize the tangent edge lengths.

## Conclusion

The v2 diagnostics separate three conditions that v1 could not distinguish:

1. intrinsic curvature closure;
2. coordinate realization of tangent edge lengths;
3. coordinate-level theta agreement.

For the production runs above:

- `run_16_unbounded_embed` reaches an intrinsic Ricci endpoint, but not a coordinate-level strong-INC candidate;
- bounded runs keep radii more reasonable but block Ricci convergence under the requested bounds;
- none of the tested randomized production runs satisfy the coordinate-level strong-INC candidate thresholds.

A valid coordinate-level strong-INC candidate still requires all three:

```text
max_abs_K < tol
max_edge_abs_residual < embed_tol
theta_mean_abs < 1e-6
```

No production run met all three criteria.
