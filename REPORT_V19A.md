# INC-Ricci v19A overlap-aware continuation report

## Build And Tests

Build under WSL with voro++:

```bash
g++ -O3 -std=c++17 inc_ricci_v19A.cpp -o inc_ricci_v19A \
  -I/home/tianyu/miniconda3/envs/analyse/include \
  -I/home/tianyu/miniconda3/envs/analyse/include/voro++ \
  -L/home/tianyu/miniconda3/envs/analyse/lib \
  -lvoro++ \
  -Wl,-rpath,/home/tianyu/miniconda3/envs/analyse/lib \
  -pthread
```

Test command:

```bash
./inc_ricci_v19A --test
```

Result: all v19A tests passed.

| Test | Result | Key diagnostic |
|---|---|---|
| Hex baseline | PASS | `Theta_sum=2.41e-16`, `Theta_shape=3.70e-16`, `du_rms=2.76e-22` |
| Shape analytic gradient | PASS | worst relative error `2.02e-10` |
| Benchmark A small | PASS | `Theta_shape 3.65e-02 -> 1.67e-03` |
| Benchmark B small | PASS | `Theta_sum 2.13e-02 -> 3.89e-03` |
| Limit-cycle diagnostic | PASS | artificial zig-zag `D_net=0` |
| F decomposition | PASS | component sum equals `F_total` |
| Pareto archive | PASS | archive produced |

## Implementation Summary

`inc_ricci_v19A.cpp` preserves v18 operator blocks and adds:

- cycle-level `F_total` decomposition into `F_sum`, `F_shape`, `F_U`, `F_radius`, and `F_invalid`;
- explicit cycle rejection reasons;
- Metropolis cycle acceptance with optional auto-temperature;
- overlap-weight continuation schedule over `w_U`;
- adaptive `bounded_du`;
- Pareto archive over theta, overlap, hard-response, and radius-cost objectives;
- optional contact-weighted shape edge target;
- `plot_v19A_results.py`.

This version is diagnostic. It does not claim INC generation unless theta defects and overlap are all low.

## N=256 Cycle Proposal Diagnostic

Command used 100 proposals at `phi=0.84`, `bounded_du=0.03`, `shape_gamma=0.05`, `w_U=0.5`.

| Statistic | Value |
|---|---:|
| hard-valid fraction | `1.0` |
| downhill fraction | `1.0` |
| uphill-but-hard-valid fraction | `0.0` |
| dominant F component | `F_shape` |

Interpretation: for these v19A settings, cycle proposals are not blocked by hard constraints. Annealing is not needed for the sampled proposals because they are already downhill in the configured objective.

## N=256 Overlap Continuation

Command used `w_U_list=0.05,0.10,0.20,0.50,1.0`, `100` cycles per stage, adaptive `bounded_du`.

| Diagnostic | Initial | Final |
|---|---:|---:|
| `Theta_sum` | `0.0665912` | `0.0622028` |
| `Theta_shape` | `0.117539` | `0.0846684` |
| `Theta_NM` | `0.117236` | `0.0972383` |
| `U_overlap` | `3.273e-04` | `0.167644` |
| `du_rms` | `0.131671` | `0.124347` |
| accepted cycles | - | `54 / 500` |

The run produced diagonal improvement in the defect plane, but the improvement is paid for by a large overlap increase.

## N=256 Contact-Weighted Edge Target

Same continuation schedule, with `--shape_edge_mode contact_weighted --shape_gamma 0.10 --shape_g0_factor 0.05`.

| Diagnostic | Initial | Final |
|---|---:|---:|
| `Theta_sum` | `0.0665912` | `0.0622028` |
| `Theta_shape` | `0.117539` | `0.0867198` |
| `Theta_NM` | `0.117236` | `0.0982286` |
| `U_overlap` | `3.273e-04` | `0.152472` |
| `du_rms` | `0.131671` | `0.124347` |
| accepted cycles | - | `49 / 500` |

Contact weighting slightly reduces the overlap cost compared with global-gamma continuation, but gives a slightly weaker shape-channel improvement.

## N=576 Scaling Check

Command used `nx=24`, `ny=24`, `radius_seed=202`, `relax_steps=30000`, and the same continuation schedule.

| Diagnostic | Initial | Final |
|---|---:|---:|
| `Theta_sum` | `0.0774266` | `0.0723808` |
| `Theta_shape` | `0.112022` | `0.0928139` |
| `Theta_NM` | `0.116785` | `0.104404` |
| `U_overlap` | `5.446e-04` | `0.134268` |
| `du_rms` | `0.151646` | `0.143061` |
| accepted cycles | - | `43 / 500` |
| runtime | - | about `1033 s` |

The N=576 run preserves the qualitative N=256 behavior: overlap-aware continuation opens a diagonal path in defect space, but does not recover low overlap by the final high-`w_U` stage.

## Pareto Interpretation

For N=256 contact-weighted continuation:

- best `Theta_NM`: `0.09823` at `U_overlap=0.15247`;
- best low-overlap point remains the initial state with `U_overlap=3.27e-04`;
- compromise with overlap weight `c=1` selects an early state: `Theta_NM=0.11073`, `U_overlap=0.00590`.

For N=576:

- best `Theta_NM`: `0.10440` at `U_overlap=0.13427`;
- best low-overlap point remains the initial state with `U_overlap=5.45e-04`;
- compromise with `c=0.1` selects an early state: `Theta_NM=0.10923`, `U_overlap=0.01796`.

## Conclusion

v19A answers the immediate v18 failure-forensics question:

```text
cycle-level freeze was protocol/weight related: yes
hard constraints were the main blocker under v19A settings: no
overlap-aware continuation opens diagonal motion: yes
high-w_U stages recover low overlap: no
```

The continuation protocol can reduce both `Theta_sum` and `Theta_shape`, and this persists at N=576. However, the final states retain substantial overlap. The current result is therefore a Pareto tradeoff, not an INC endpoint.

The next useful direction is to add a mechanism that relaxes overlap after geometric ordering without undoing the theta gains, or to make overlap part of a staged compression/decompression protocol rather than a static penalty.
