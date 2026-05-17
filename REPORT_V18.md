# INC-Ricci v18 operator-validation report

## Build And Tests

Build command used under WSL:

```bash
g++ -O3 -std=c++17 inc_ricci_v18.cpp -o inc_ricci_v18 \
  -I/home/tianyu/miniconda3/envs/analyse/include \
  -I/home/tianyu/miniconda3/envs/analyse/include/voro++ \
  -L/home/tianyu/miniconda3/envs/analyse/lib \
  -lvoro++ \
  -Wl,-rpath,/home/tianyu/miniconda3/envs/analyse/lib \
  -pthread
```

Test command:

```bash
./inc_ricci_v18 --test
```

Result: all v18 tests passed.

| Test | Result | Key diagnostics |
|---|---|---|
| Hex baseline | PASS | `Theta_sum=2.41e-16`, `Theta_shape=3.70e-16`, `du_rms=2.76e-22` |
| Shape analytic gradient | PASS | worst relative error `2.02e-10` |
| Benchmark A small | PASS | `Theta_shape 3.65e-02 -> 1.67e-03` |
| Benchmark B small | PASS | `Theta_sum 2.13e-02 -> 3.89e-03` |
| Limit-cycle diagnostic | PASS | artificial zig-zag gives `D_net=0` |

## v18 Implementation Summary

`inc_ricci_v18.cpp` preserves the v17/v16 voro++ radical topology and defect-plane diagnostics, and adds validated operator blocks:

- analytic-gradient coordinate shape block using cosine angle energy;
- metric-coupled edge regularization with tunable `shape_gamma`;
- bounded / regularized Ricci block;
- finite-difference box aspect block;
- topology surgery with micro-quench;
- coupled cycle protocol with block-level response matrix;
- limit-cycle diagnostic;
- controlled benchmarks A/B/C;
- polydisperse application mode.

The v18 state is treated as:

```text
S = (X, u, T, H)
```

where `X` is coordinate state, `u=log(r)`, `T` is topology, and `H` is the periodic box geometry.

## Benchmark A: Shape Operator Validation

Command:

```bash
./inc_ricci_v18 --nx 16 --ny 16 \
  --box_mode rectangular_hex \
  --benchmark_A_coord_noise_hex 1 \
  --coord_noise_amp 0.05 \
  --shape_gamma_scan 0,0.1,0.25,0.5,0.75,1.0 \
  --shape_steps 1000 \
  --out v18_A_gamma_scan \
  --progress_every 0
```

| `shape_gamma` | `Theta_shape` initial | `Theta_shape` final | `Theta_sum` final | `U_final` | success |
|---:|---:|---:|---:|---:|---:|
| 0.00 | `3.4935e-02` | `1.6619e-02` | `2.55e-16` | `2.2252e-01` | 1 |
| 0.10 | `3.4935e-02` | `1.4950e-02` | `2.57e-16` | `1.8028e-01` | 1 |
| 0.25 | `3.4935e-02` | `1.2448e-02` | `2.63e-16` | `1.2523e-01` | 1 |
| 0.50 | `3.4935e-02` | `8.2869e-03` | `2.55e-16` | `5.5682e-02` | 1 |
| 0.75 | `3.4935e-02` | `4.1377e-03` | `2.50e-16` | `1.3926e-02` | 1 |
| 1.00 | `3.4935e-02` | `1.3232e-13` | `2.92e-16` | `7.77e-24` | 1 |

Best gamma: `1.0`.

Interpretation: the analytic shape operator is mathematically valid on the coordinate-noise hex benchmark. Full metric coupling (`gamma=1`) repairs coordinate angular disorder best while keeping `Theta_sum` and `du_rms` at the zero baseline.

## Benchmark B: Bounded Ricci Validation

Command:

```bash
./inc_ricci_v18 --nx 16 --ny 16 \
  --box_mode rectangular_hex \
  --benchmark_B_radius_noise_hex 1 \
  --radius_noise_amp 0.05 \
  --bounded_du 0.10 \
  --ricci_steps 500 \
  --ricci_dt 0.01 \
  --out v18_B_ricci_hex \
  --progress_every 0
```

| Diagnostic | Initial | Final |
|---|---:|---:|
| `Theta_sum` | `2.4615e-02` | `8.8635e-03` |
| `Theta_shape` | `1.3783e-02` | `7.3576e-03` |
| `du_rms` | `5.1642e-02` | `2.6141e-02` |
| `max |u-u0|` | - | `1.0000e-01` |

Interpretation: the bounded Ricci block reduces the radius-topology channel and the hard proxy while respecting the configured trust-region boundary. It reaches the `bounded_du=0.10` limit.

## Benchmark C: Coupled Repair

Command:

```bash
./inc_ricci_v18 --nx 16 --ny 16 \
  --box_mode rectangular_hex \
  --benchmark_C_mixed_noise_hex 1 \
  --coord_noise_amp 0.05 \
  --radius_noise_amp 0.05 \
  --cycles 200 \
  --shape_steps_per_cycle 50 \
  --ricci_steps_per_cycle 50 \
  --shape_gamma 1.0 \
  --bounded_du 0.10 \
  --out v18_C_coupled_hex \
  --progress_every 0
```

| Diagnostic | Initial | Final |
|---|---:|---:|
| `Theta_sum` | `2.4615e-02` | `8.8635e-03` |
| `Theta_shape` | `3.7456e-02` | `1.2926e-02` |
| `Theta_NM` | `4.4594e-02` | `1.5593e-02` |
| `B_rms` | `1.8220e-01` | `6.5156e-02` |
| `du_rms` | `5.1642e-02` | `2.6141e-02` |
| `U_overlap` | `2.2229` | `1.8406e-01` |
| `D_net` | - | `0.7545` |
| zigzag index | - | `0.2455` |

Interpretation: the controlled mixed-noise case repairs diagonally in the defect plane. The accepted block trajectory is directed enough (`D_net=0.7545`) that the benchmark is not dominated by limit-cycle behavior.

## Polydisperse Application

Command:

```bash
./inc_ricci_v18 --nx 16 --ny 16 \
  --radius_seed 101 \
  --target_dist lognormal \
  --target_delta 0.25 \
  --target_radius_ratio 5 \
  --coord_mode soft_relax \
  --target_phi 0.84 \
  --relax_steps 20000 \
  --use_voropp 1 \
  --compute_nm_theta 1 \
  --decompose_theta 1 \
  --apply_to_polydisperse 1 \
  --cycles 200 \
  --shape_steps_per_cycle 50 \
  --ricci_steps_per_cycle 50 \
  --shape_gamma 1.0 \
  --bounded_du 0.10 \
  --out v18_polydisperse_16 \
  --progress_every 50
```

| Diagnostic | Initial | Final |
|---|---:|---:|
| `Theta_sum` | `6.6591e-02` | `6.6591e-02` |
| `Theta_shape` | `1.1754e-01` | `1.1754e-01` |
| `Theta_NM` | `1.1724e-01` | `1.1724e-01` |
| `B_rms` | `4.8845e-01` | `4.8845e-01` |
| `du_rms` | `1.3167e-01` | `1.3167e-01` |
| `U_overlap` | `3.2731e-04` | `3.2731e-04` |
| `D_net` | - | `nan` |

Interpretation: the validated blocks do not accept a net improvement on the real N=256 soft-relaxed polydisperse radical sample under the current trust-region/cycle protocol. The accepted path length is zero, so `D_net` is `nan`.

## Conclusion

v18 validates the operator blocks on controlled baselines:

```text
shape block works on coordinate-noise hex: yes
bounded Ricci works on radius-noise hex: yes
coupled cycle repairs mixed-noise hex diagonally: yes
same protocol improves real polydisperse sample: no
```

Current interpretation:

- The analytic coordinate-gradient shape operator is correct and effective on a clean shape-only defect.
- The bounded Ricci operator is effective on a clean radius-topology defect.
- The coupled protocol works when the two channels are synthetically separable.
- The real polydisperse radical topology remains frustrated under the simplified v18 trust-region dynamics.

This means the immediate bottleneck is no longer basic operator correctness. The missing ingredient is likely the coupling mechanism that creates accepted topology/coordinate/radius moves in genuinely amorphous polydisperse states.
